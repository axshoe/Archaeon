"""
src/analysis/scorer.py
----------------------
Industrial Applicability Score (IAS) computation for Archaeon.

The IAS is a composite score (0-100) that ranks candidate enzyme sequences
by their predicted value for industrial applications. It combines three
independent signal categories:

  1. Thermostability Features (40% weight)
     Sequence-based thermostability signals from features.py.
     Includes aliphatic index, instability index, charged residue ratio,
     proline content, and aromaticity. These are the features with
     the strongest literature support for thermostability prediction.

  2. Sequence Quality and Length (20% weight)
     Sequence length in the industrial "sweet spot" and absence of
     ambiguous residues. Very short sequences (<150 aa) are likely
     fragments; very long sequences (>1500 aa) are less tractable as
     standalone industrial enzymes.

  3. BLAST Identity to Reference Enzymes (40% weight)
     Percent identity to a known well-characterized thermostable enzyme
     from the same family. Higher identity = more confident the candidate
     is a functional variant of a known industrial enzyme.
     If BLAST was not run (to save time), this component defaults to 50/100.

Scoring philosophy:
  This score is not a thermostability prediction; it is an industrial
  applicability ranking. It is designed to surface candidates that are
  (a) likely thermostable, (b) likely to be real enzymes with known function,
  and (c) in a tractable size range for recombinant expression.

  The weights (40/20/40) reflect practical industrial priorities. BLAST
  identity gets 40% because it provides function annotation, which is the
  single most important factor for industrial enzyme development. You cannot
  develop an enzyme if you don't know what it does.

  The scores are NOT calibrated against experimental data; we don't have
  a training set. They are theoretically grounded and designed for relative
  ranking within the candidate set, not for absolute prediction.
"""

import math
import logging
from typing import Optional

from src.analysis.features import extract_all_features, thermostability_summary

logger = logging.getLogger(__name__)

# -------------------------------------------------------------------------
# IAS Component Weights
# Must sum to 1.0
# -------------------------------------------------------------------------
W_THERMOSTABILITY = 0.40  # Sequence thermostability features
W_QUALITY         = 0.20  # Sequence quality and length
W_BLAST           = 0.40  # BLAST identity to reference enzyme


def score_thermostability_component(features: dict) -> float:
    """
    Compute the thermostability component score (0-100) from sequence features.

    Uses a weighted combination of features with literature-informed weights:
      - Aliphatic index: 30% (strongest single predictor of thermostability)
      - Instability index: 25% (inversely related; lower II = more stable)
      - Total charged fraction: 25% (ion pair / salt bridge density signal)
      - Proline content: 10% (loop rigidity signal)
      - Aromaticity: 10% (aromatic stacking network signal)

    Each feature is normalized to a 0-1 scale based on biological reference ranges
    derived from the extremophile enzyme literature. The normalization ranges
    are documented inline.

    Returns a score in [0, 100].
    """
    if not features:
        return 0.0

    # --- Aliphatic index normalization ---
    # Reference range from thermostable enzyme literature:
    #   < 60: mesophilic (score = 0)
    #   60-80: moderate thermostability
    #   > 120: highly thermostable (score = 1.0)
    ai = features.get("aliphatic_index", 0)
    ai_norm = min(max((ai - 60) / 60.0, 0.0), 1.0)  # maps 60-120 to 0-1

    # --- Instability index normalization ---
    # INVERTED: lower instability index = more stable = higher score
    #   > 70: very unstable (score = 0)
    #   < 20: very stable (score = 1.0)
    # Typical thermostable enzymes: II between 20-40
    ii = features.get("instability_index", 70.0)
    ii_norm = min(max((70.0 - ii) / 50.0, 0.0), 1.0)  # maps 20-70 to 1-0 inverted

    # --- Total charged fraction normalization ---
    # Reference range:
    #   < 0.10: low (mesophilic typical)
    #   0.20-0.35: thermophile optimal range
    #   > 0.40: extreme (may affect solubility)
    tcf = features.get("total_charged_fraction", 0.0)
    # Penalize both extremes; peak at 0.25
    tcf_norm = 1.0 - min(abs(tcf - 0.25) / 0.25, 1.0)

    # --- Proline content normalization ---
    # Thermostable enzymes: 0.04-0.10 proline fraction
    # > 0.10 can indicate structural issues
    pro = features.get("proline_fraction", 0.0)
    pro_norm = min(max((pro - 0.02) / 0.08, 0.0), 1.0)

    # --- Aromaticity normalization ---
    # Typical range: 0.05-0.15 for thermostable proteins
    aro = features.get("aromaticity", 0.0)
    aro_norm = min(max((aro - 0.05) / 0.10, 0.0), 1.0)

    # Weighted sum of component scores
    thermo_score = (
        0.30 * ai_norm +
        0.25 * ii_norm +
        0.25 * tcf_norm +
        0.10 * pro_norm +
        0.10 * aro_norm
    ) * 100.0

    return round(thermo_score, 2)


def score_quality_component(features: dict) -> float:
    """
    Compute sequence quality and length score (0-100).

    Length penalty curve:
      - Length in [150, 800]: score 100 (industrial sweet spot)
      - Length in [100, 150) or (800, 1200]: score decays linearly
      - Length < 100 or > 1200: score 0

    Rationale: Industrial enzyme applications prefer proteins in the 150-800 aa
    range. Shorter sequences are likely not full-length enzymes; they may be
    catalytic domains that lack the accessory domains needed for proper folding.
    Longer sequences (>800 aa) are harder to express in E. coli (a workhorse
    production host) and may be multifunctional complexes.
    """
    if not features:
        return 0.0

    length = features.get("length", 0)

    # Length score using piecewise linear function
    if 150 <= length <= 800:
        length_score = 1.0
    elif 100 <= length < 150:
        length_score = (length - 100) / 50.0
    elif 800 < length <= 1200:
        length_score = 1.0 - (length - 800) / 400.0
    else:
        length_score = 0.0

    return round(length_score * 100.0, 2)


def score_blast_component(blast_hits: Optional[list]) -> tuple:
    """
    Compute BLAST identity component score (0-100).

    Args:
        blast_hits: List of hit dicts from ncbi.blast_against_enzyme_families(),
                    or None if BLAST was skipped.

    Returns:
        (blast_score, best_family, best_identity_pct)

    Score mapping:
      - Identity >= 80%: score 100 (very close homolog, confident function annotation)
      - Identity 60-80%: score scales 70-100 (likely same enzyme family)
      - Identity 40-60%: score scales 30-70 (distant homolog, uncertain function)
      - Identity < 40%: score < 30 (very remote homology)
      - No BLAST data: score 50 (neutral, no evidence either way)

    This is the most important component because function annotation is the
    primary industrial bottleneck. A candidate with 85% identity to a known
    thermostable cellulase from Thermotoga maritima is immediately actionable
    for industrial cellulose degradation. A candidate with 25% identity to
    anything is a much larger experimental investment.
    """
    if not blast_hits:
        # BLAST skipped or failed; return neutral score with no annotation
        return 50.0, "unknown", 0.0

    # Take the best hit (already sorted by identity in ncbi.py)
    best_hit = blast_hits[0]
    identity = best_hit.get("identity_pct", 0.0)
    family = best_hit.get("family", "unknown")

    # Piecewise linear identity -> score mapping
    if identity >= 80:
        blast_score = 100.0
    elif identity >= 60:
        blast_score = 70.0 + (identity - 60) / 20.0 * 30.0
    elif identity >= 40:
        blast_score = 30.0 + (identity - 40) / 20.0 * 40.0
    elif identity >= 20:
        blast_score = 10.0 + (identity - 20) / 20.0 * 20.0
    else:
        blast_score = identity / 20.0 * 10.0

    return round(blast_score, 2), family, identity


def compute_ias(
    sequence: str,
    blast_hits: Optional[list] = None,
    organism: str = "unknown",
    source_biome: str = "unknown",
) -> dict:
    """
    Compute the full Industrial Applicability Score for a single candidate.

    This is the main entry point for scoring. Returns a comprehensive dict
    with all intermediate scores and the final IAS, suitable for direct
    inclusion in the report and CSV output.

    Args:
        sequence:     Amino acid sequence string
        blast_hits:   BLAST results from ncbi.blast_against_enzyme_families()
                      (None if BLAST was not run)
        organism:     Source organism name (metadata only; not used in scoring)
        source_biome: Isolation environment (metadata only)

    Returns:
        {
          "ias": float,                     # Final composite score [0, 100]
          "thermostability_score": float,   # Component 1 score [0, 100]
          "quality_score": float,           # Component 2 score [0, 100]
          "blast_score": float,             # Component 3 score [0, 100]
          "predicted_family": str,          # Predicted enzyme family
          "best_blast_identity": float,     # Best BLAST % identity
          "length": int,                    # Sequence length
          "features": dict,                 # All extracted sequence features
          "thermostability_signals": list,  # Human-readable feature interpretation
          "organism": str,
          "source_biome": str,
        }
    """
    features = extract_all_features(sequence)
    if not features:
        logger.warning("Could not extract features from sequence; returning zero score.")
        return {
            "ias": 0.0,
            "thermostability_score": 0.0,
            "quality_score": 0.0,
            "blast_score": 0.0,
            "predicted_family": "unknown",
            "best_blast_identity": 0.0,
            "length": len(sequence),
            "features": {},
            "thermostability_signals": [],
            "organism": organism,
            "source_biome": source_biome,
        }

    # Compute three component scores
    thermo_score  = score_thermostability_component(features)
    quality_score = score_quality_component(features)
    blast_score, predicted_family, best_identity = score_blast_component(blast_hits)

    # Weighted composite: IAS = 0.40 * thermo + 0.20 * quality + 0.40 * blast
    ias = (
        W_THERMOSTABILITY * thermo_score +
        W_QUALITY         * quality_score +
        W_BLAST           * blast_score
    )

    # Get qualitative thermostability signals for the report narrative
    ts_summary = thermostability_summary(sequence)
    thermo_signals = ts_summary.get("signals", [])

    return {
        "ias":                    round(ias, 2),
        "thermostability_score":  round(thermo_score, 2),
        "quality_score":          round(quality_score, 2),
        "blast_score":            round(blast_score, 2),
        "predicted_family":       predicted_family,
        "best_blast_identity":    best_identity,
        "length":                 features.get("length", 0),
        "features":               features,
        "thermostability_signals": thermo_signals,
        "organism":               organism,
        "source_biome":           source_biome,
    }


def rank_candidates(
    records: list,
    blast_results: Optional[dict] = None,
    top_n: int = 20,
) -> list:
    """
    Score and rank a list of sequence records by IAS.

    Args:
        records:      List of dicts from ncbi.py or mgnify.py, each having
                      at least 'sequence', 'organism', 'source_biome', 'id'.
        blast_results: Dict mapping record id -> blast_hits list. None if BLAST
                       was not run.
        top_n:        Return only the top N candidates after ranking.

    Returns:
        Sorted list of result dicts (highest IAS first), each record's original
        metadata merged with the IAS scoring dict.
    """
    if not records:
        logger.warning("No records to score.")
        return []

    logger.info(f"Scoring {len(records)} candidates...")
    scored = []

    for rec in records:
        seq = rec.get("sequence", "")
        if not seq:
            continue

        blast_hits = None
        if blast_results:
            blast_hits = blast_results.get(rec.get("id", ""), None)

        ias_result = compute_ias(
            sequence=seq,
            blast_hits=blast_hits,
            organism=rec.get("organism", "unknown"),
            source_biome=rec.get("source_biome", "unknown"),
        )

        # Merge record metadata with scoring result
        full_result = {**rec, **ias_result}
        scored.append(full_result)

    # Sort by IAS descending
    scored.sort(key=lambda x: x.get("ias", 0), reverse=True)

    top = scored[:top_n]
    logger.info(
        f"Ranking complete. Top IAS: {top[0]['ias']:.1f} | "
        f"Lowest in top-{top_n}: {top[-1]['ias']:.1f}"
    )
    return top
