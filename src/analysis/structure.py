"""
src/analysis/structure.py
-------------------------
ESMFold API integration for 3D protein structure prediction.

ESMFold is Meta AI's protein structure prediction model (Lin et al., 2022,
Science 379:1123-1130). Unlike AlphaFold2, ESMFold:
  - Does NOT require a multiple sequence alignment (MSA)
  - Runs in seconds per sequence (vs. minutes for AlphaFold2 with MSA)
  - Is available via a free public REST API at ESM Metagenomic Atlas

This makes ESMFold ideal for metagenomic sequences, which often have few
close homologs in sequence databases (making MSA construction poor anyway).
The tradeoff: ESMFold is somewhat less accurate than AlphaFold2 on sequences
with many known homologs. For novel extremophile sequences from undersampled
biomes, this tradeoff is acceptable.

API endpoint: https://api.esmatlas.com/foldSequence/v1/pdb/
Returns: PDB format structure file as plain text.

pLDDT score: ESMFold outputs a per-residue confidence score (pLDDT, 0-100)
embedded in the B-factor column of the PDB file. High pLDDT (>70) indicates
confident structural prediction; low pLDDT (<50) indicates disordered or
poorly predicted regions.

This module handles:
  1. Submitting sequences to ESMFold API
  2. Parsing the PDB response to extract pLDDT scores
  3. Computing mean pLDDT and per-region quality summaries
  4. Saving PDB files for visualization in py3Dmol
"""

import time
import logging
import re
from pathlib import Path
from typing import Optional

import requests

logger = logging.getLogger(__name__)

ESMFOLD_API_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"

# ESMFold API limits
ESMFOLD_MAX_LENGTH = 400   # Sequences >400 aa will be truncated with a warning
ESMFOLD_REQUEST_INTERVAL = 2.0  # Seconds between requests (be polite to EBI)

# pLDDT quality thresholds (same scale as AlphaFold2)
PLDDT_VERY_HIGH = 90   # Very high confidence
PLDDT_HIGH      = 70   # Confident prediction
PLDDT_LOW       = 50   # Low confidence; likely disordered
PLDDT_VERY_LOW  = 30   # Very low confidence; do not interpret


def predict_structure(
    sequence: str,
    candidate_id: str,
    output_dir: Optional[str] = None,
    truncate: bool = True,
) -> dict:
    """
    Submit a protein sequence to ESMFold API and retrieve the predicted structure.

    Args:
        sequence:     Amino acid sequence (standard 1-letter code)
        candidate_id: Unique identifier for this candidate (used for filename)
        output_dir:   Directory to save the PDB file; if None, PDB not saved
        truncate:     If True, truncate sequences >ESMFOLD_MAX_LENGTH aa with warning.
                      If False, raise ValueError on overlong sequences.

    Returns:
        {
          "candidate_id": str,
          "pdb_string": str,         # Full PDB file as string
          "pdb_path": str or None,   # Path where PDB was saved
          "mean_plddt": float,       # Mean per-residue confidence (0-100)
          "plddt_scores": list,      # Per-residue pLDDT values
          "quality_summary": dict,   # Fraction of residues in each pLDDT tier
          "sequence_length": int,    # Length of sequence actually predicted
          "success": bool,
          "error": str or None,
        }
    """
    result = {
        "candidate_id":  candidate_id,
        "pdb_string":    None,
        "pdb_path":      None,
        "mean_plddt":    0.0,
        "plddt_scores":  [],
        "quality_summary": {},
        "sequence_length": len(sequence),
        "success":       False,
        "error":         None,
    }

    # Length check
    if len(sequence) > ESMFOLD_MAX_LENGTH:
        if truncate:
            logger.warning(
                f"Sequence {candidate_id} has length {len(sequence)} > "
                f"{ESMFOLD_MAX_LENGTH}. Truncating to N-terminal {ESMFOLD_MAX_LENGTH} aa. "
                "This affects pLDDT interpretation; only the truncated region is predicted."
            )
            sequence = sequence[:ESMFOLD_MAX_LENGTH]
            result["sequence_length"] = ESMFOLD_MAX_LENGTH
        else:
            result["error"] = f"Sequence too long: {len(sequence)} aa (max {ESMFOLD_MAX_LENGTH})"
            return result

    # Validate sequence characters (only standard amino acids)
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
    invalid = set(sequence.upper()) - valid_aa
    if invalid:
        logger.warning(f"Non-standard residues in {candidate_id}: {invalid}. Filtering.")
        sequence = "".join(aa for aa in sequence.upper() if aa in valid_aa)
        if len(sequence) < 10:
            result["error"] = "Too few valid residues after filtering"
            return result

    logger.info(f"Submitting {candidate_id} (length {len(sequence)}) to ESMFold API")

    # Submit to ESMFold API
    # The API accepts raw sequence as plain text in the POST body.
    try:
        response = requests.post(
            ESMFOLD_API_URL,
            data=sequence,
            headers={"Content-Type": "application/x-www-form-urlencoded"},
            timeout=120,  # ESMFold can be slow for longer sequences
        )
        if response.status_code != 200:
            result["error"] = f"ESMFold API returned HTTP {response.status_code}: {response.text[:200]}"
            logger.error(result["error"])
            return result

        pdb_string = response.text

    except requests.RequestException as e:
        result["error"] = f"ESMFold API request failed: {e}"
        logger.error(result["error"])
        return result

    # Parse pLDDT scores from PDB B-factor column.
    # In ESMFold PDB output, the B-factor (column 61-66 in ATOM records)
    # contains the per-residue pLDDT score. We parse one score per residue
    # (take the CA atom to avoid duplicating per residue).
    plddt_scores = _parse_plddt_from_pdb(pdb_string)

    if not plddt_scores:
        result["error"] = "Could not parse pLDDT scores from PDB output"
        logger.error(result["error"])
        return result

    mean_plddt = sum(plddt_scores) / len(plddt_scores)
    quality_summary = _compute_quality_summary(plddt_scores)

    # Save PDB file if output directory specified
    pdb_path = None
    if output_dir:
        pdb_dir = Path(output_dir)
        pdb_dir.mkdir(parents=True, exist_ok=True)
        pdb_path = str(pdb_dir / f"{candidate_id}.pdb")
        with open(pdb_path, "w") as f:
            f.write(pdb_string)
        logger.info(f"Saved PDB to {pdb_path}")

    result.update({
        "pdb_string":      pdb_string,
        "pdb_path":        pdb_path,
        "mean_plddt":      round(mean_plddt, 2),
        "plddt_scores":    plddt_scores,
        "quality_summary": quality_summary,
        "success":         True,
        "error":           None,
    })

    time.sleep(ESMFOLD_REQUEST_INTERVAL)  # Rate limit: be polite to EBI servers
    return result


def _parse_plddt_from_pdb(pdb_string: str) -> list:
    """
    Extract per-residue pLDDT scores from ESMFold PDB output.

    ESMFold stores pLDDT in the B-factor column of ATOM records.
    We use only CA (alpha-carbon) atoms to get one value per residue,
    avoiding duplicate entries for side-chain atoms.

    PDB ATOM record format (fixed-width):
      Columns 1-6:   Record type ("ATOM  ")
      Columns 13-16: Atom name (e.g., " CA ")
      Columns 61-66: B-factor (pLDDT here)

    Returns list of float pLDDT values, one per residue.
    """
    plddt_scores = []
    for line in pdb_string.split("\n"):
        if not line.startswith("ATOM"):
            continue
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue
        try:
            b_factor = float(line[60:66].strip())
            if b_factor <= 1.0:
                b_factor *= 100.0
            plddt_scores.append(b_factor)
        except (ValueError, IndexError):
            continue

    return plddt_scores


def _compute_quality_summary(plddt_scores: list) -> dict:
    """
    Compute the fraction of residues in each pLDDT quality tier.

    Returns:
      {
        "very_high": float,  # Fraction with pLDDT >= 90
        "high": float,       # Fraction with 70 <= pLDDT < 90
        "low": float,        # Fraction with 50 <= pLDDT < 70
        "very_low": float,   # Fraction with pLDDT < 50
        "n_residues": int,
      }
    """
    if not plddt_scores:
        return {}

    n = len(plddt_scores)
    very_high = sum(1 for s in plddt_scores if s >= PLDDT_VERY_HIGH) / n
    high      = sum(1 for s in plddt_scores if PLDDT_HIGH <= s < PLDDT_VERY_HIGH) / n
    low       = sum(1 for s in plddt_scores if PLDDT_LOW  <= s < PLDDT_HIGH) / n
    very_low  = sum(1 for s in plddt_scores if s < PLDDT_LOW) / n

    return {
        "very_high":  round(very_high, 3),
        "high":       round(high, 3),
        "low":        round(low, 3),
        "very_low":   round(very_low, 3),
        "n_residues": n,
    }


def batch_predict_structures(
    candidates: list,
    output_dir: str,
    max_candidates: int = 10,
) -> dict:
    """
    Run ESMFold structure prediction on multiple candidates.

    Only processes the top `max_candidates` sequences to keep runtime
    manageable. At ~5 seconds per sequence, 10 candidates = ~50 seconds.

    Args:
        candidates:     List of scoring dicts from scorer.rank_candidates(),
                        each must have 'id' and 'sequence' keys.
        output_dir:     Directory for PDB file output
        max_candidates: Maximum number of structure predictions to run

    Returns:
        Dict mapping candidate_id -> structure prediction result dict
    """
    results = {}
    target = candidates[:max_candidates]

    logger.info(f"Running ESMFold on top {len(target)} candidates")

    for i, candidate in enumerate(target):
        cid = candidate.get("id", f"candidate_{i}")
        seq = candidate.get("sequence", "")

        if not seq:
            logger.warning(f"No sequence for candidate {cid}; skipping")
            continue

        logger.info(f"  [{i+1}/{len(target)}] Predicting structure for {cid}")
        result = predict_structure(
            sequence=seq,
            candidate_id=cid,
            output_dir=output_dir,
        )
        results[cid] = result

        if result["success"]:
            logger.info(
                f"    Mean pLDDT: {result['mean_plddt']:.1f} | "
                f"High confidence: {result['quality_summary'].get('very_high', 0)*100:.0f}% of residues"
            )
        else:
            logger.warning(f"    Failed: {result['error']}")

    successes = sum(1 for r in results.values() if r["success"])
    logger.info(f"Structure prediction complete: {successes}/{len(target)} succeeded")
    return results


def plddt_to_color_hex(plddt: float) -> str:
    """
    Convert a pLDDT score to a color hex code for visualization.
    Matches AlphaFold2/ESMFold standard color scheme:
      >= 90: #0053D6 (dark blue, very high confidence)
      70-90: #65CBF3 (light blue, confident)
      50-70: #FFDB13 (yellow, low confidence)
      < 50:  #FF7D45 (orange, very low confidence)
    """
    if plddt >= 90:
        return "#0053D6"
    elif plddt >= 70:
        return "#65CBF3"
    elif plddt >= 50:
        return "#FFDB13"
    else:
        return "#FF7D45"
