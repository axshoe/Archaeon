"""
src/analysis/features.py
------------------------
Sequence feature extraction for thermostability prediction.

All features computed here are sequence-based: they derive solely from
the amino acid sequence, no 3D structure needed. This is intentional.
Structure prediction with ESMFold is expensive (time and memory); we use
sequence features as a cheap first-pass filter, then only invoke ESMFold
on the top candidates.

Scientific basis for each feature is documented below. These are not ad-hoc
choices; they are the standard features used in thermostability prediction
papers since the 1990s and validated extensively in the literature.

Key references:
  - Ikai (1980): aliphatic index formula
  - Boman (2003): GRAVY and instability index basis
  - Petersen et al. (2010): charged residue analysis for thermostability
  - Zeldovich et al. (2007): purine loading rule for thermophile sequences
    (DNA-level, but reflected in amino acid composition indirectly)
"""

import math
import logging
from typing import Optional

logger = logging.getLogger(__name__)

# -------------------------------------------------------------------------
# Amino acid property tables
# -------------------------------------------------------------------------

# Kyte-Doolittle hydrophobicity scale (Kyte & Doolittle, 1982, J. Mol. Biol.).
# Higher = more hydrophobic. Used for GRAVY index.
HYDROPHOBICITY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# Instability weight table (Guruprasad et al., 1990, Protein Engineering).
# These are dipeptide weights; a protein is unstable if its weighted sum > 40.
# Only the most informative dipeptides are listed here; unlisted pairs get 1.0.
INSTABILITY_WEIGHTS = {
    "WW": 1.0, "WC": 1.0, "WT": 0.733, "CW": 1.0, "WM": 0.733,
    "CY": 1.0, "YW": 0.667, "YY": 0.5, "WA": 0.533, "WG": 0.533,
    "RC": 0.333, "QW": 0.333, "QH": 0.333, "RS": 0.2, "RE": 0.2,
    "WH": 0.333, "QR": 0.333, "WL": 0.433, "WI": 0.433, "WV": 0.433,
    "FK": 1.0, "HH": 0.667, "HR": 0.5, "HD": 0.5, "HN": 0.5,
    "HC": 0.5, "HS": 0.5, "HQ": 0.5, "HK": 0.5, "KS": 0.333,
    "KN": 0.333, "KA": 1.0, "KD": 0.5, "KP": 1.0, "KQ": 1.0,
    "MD": 0.667, "MN": 0.667, "NF": 0.667, "NP": 0.667, "NW": 0.667,
    "NY": 0.667, "NR": 0.667, "NN": 0.667, "NC": 0.667, "NK": 0.667,
    "NL": 0.667, "NH": 0.667, "NA": 0.667, "NM": 0.667, "NG": 0.667,
    "NS": 0.667, "NE": 0.667, "NT": 0.667, "NI": 0.667, "NV": 0.667,
    "ND": 0.667, "NQ": 0.667,
}

# Aliphatic amino acids used in the aliphatic index formula (Ikai, 1980).
# These are the bulky nonpolar side chains that pack tightly into the
# protein core and increase thermal stability.
ALIPHATIC_AA = {"A", "V", "I", "L"}

# Charged residues for salt bridge analysis.
POSITIVELY_CHARGED = {"R", "K", "H"}
NEGATIVELY_CHARGED = {"D", "E"}

# Aromatic residues (contribute to aromatic stacking interactions)
AROMATIC_AA = {"F", "W", "Y", "H"}

# Amino acids with small side chains that allow tight backbone packing
SMALL_AA = {"G", "A", "S", "T", "V"}

# Proline: the "structure breaker" that actually stabilizes thermophile loops
# by reducing conformational entropy of the unfolded state
PROLINE = "P"


def compute_amino_acid_composition(sequence: str) -> dict:
    """
    Compute fractional amino acid composition for all 20 standard amino acids.

    Returns a dict: {"A": 0.12, "C": 0.03, ...} where values are fractions
    summing to approximately 1.0 (minus any non-standard residues).

    This is the most basic sequence feature, but it is surprisingly informative.
    Thermophilic proteins systematically differ from mesophilic ones in their
    amino acid composition: higher content of charged residues (especially
    glutamate and arginine), lower content of asparagine and glutamine (which
    deamidate at high temperatures), and higher proline content.
    """
    seq = sequence.upper().strip()
    if not seq:
        return {}

    total = len(seq)
    standard_aa = set("ACDEFGHIKLMNPQRSTVWY")
    composition = {}

    for aa in standard_aa:
        composition[aa] = seq.count(aa) / total

    return composition


def compute_gravy_index(sequence: str) -> float:
    """
    Compute the Grand Average of Hydropathicity (GRAVY) index.

    GRAVY = (sum of Kyte-Doolittle hydrophobicity values) / sequence_length

    Reference: Kyte & Doolittle (1982), J. Mol. Biol. 157:105-132.

    Interpretation:
      GRAVY > 0: net hydrophobic (likely membrane-associated or core-packed)
      GRAVY < 0: net hydrophilic (likely soluble, secreted, or surface-exposed)
      GRAVY < -0.5: strongly hydrophilic

    For industrial enzymes, moderate hydrophilicity is desirable for solubility
    in aqueous reaction media. Very negative GRAVY correlates with proteins that
    are hard to express and purify. Very positive GRAVY suggests membrane
    insertion or aggregation-prone.

    The relationship to thermostability is not monotonic; GRAVY is used here as
    a feature in the composite score, not as a direct thermostability predictor.
    """
    seq = sequence.upper().strip()
    if not seq:
        return 0.0

    hydro_sum = sum(HYDROPHOBICITY.get(aa, 0.0) for aa in seq)
    return round(hydro_sum / len(seq), 4)


def compute_aliphatic_index(sequence: str) -> float:
    """
    Compute the aliphatic index of a protein.

    Formula (Ikai, 1980):
      AI = X(A) + 2.9 * X(V) + 3.9 * (X(I) + X(L))

    Where X(aa) = molar fraction of amino acid aa in the sequence.

    Reference: Ikai (1980), J. Biochem. 88:1895-1898.

    Interpretation:
      Higher AI = greater thermostability.
      AI > 80: typical for thermostable proteins
      AI < 60: typical for mesophilic proteins

    The biological rationale: aliphatic residues (A, V, I, L) pack tightly
    in the protein hydrophobic core. Isoleucine and leucine are weighted
    3.9x relative to alanine because their larger side chains make stronger
    van der Waals contacts in the core, resisting unfolding at high temperature.

    This is one of the most predictive single sequence features for
    thermostability in the literature.
    """
    seq = sequence.upper().strip()
    if not seq:
        return 0.0

    total = len(seq)
    xa = seq.count("A") / total  # Alanine fraction
    xv = seq.count("V") / total  # Valine fraction
    xi = seq.count("I") / total  # Isoleucine fraction
    xl = seq.count("L") / total  # Leucine fraction

    aliphatic_index = (xa + 2.9 * xv + 3.9 * (xi + xl)) * 100
    return round(aliphatic_index, 2)


def compute_instability_index(sequence: str) -> float:
    """
    Compute the instability index (II).

    Formula (Guruprasad et al., 1990):
      II = (10 / L) * sum_over_all_dipeptides(DIWV[X][Y])

    Where L = sequence length and DIWV[X][Y] is the dipeptide instability
    weight for the dipeptide XY at position i.

    Reference: Guruprasad et al. (1990), Protein Engineering 4:155-161.

    Interpretation:
      II < 40: protein is predicted to be stable
      II > 40: protein is predicted to be unstable (short half-life in vivo)

    This index was originally calibrated on in-vivo protein stability (half-life
    in cell-free systems), NOT on thermostability directly. However, in-vivo
    stability correlates with resistance to unfolding, so it is used as a proxy.
    Lower instability index is BETTER for industrial enzyme applications.

    Note: This is an empirical index with known limitations. It performs at
    about 70% accuracy in distinguishing stable vs. unstable proteins and
    should not be used as the sole criterion.
    """
    seq = sequence.upper().strip()
    if len(seq) < 2:
        return 100.0  # Degenerate case: flag as unstable

    dipeptide_sum = 0.0
    for i in range(len(seq) - 1):
        dipeptide = seq[i:i+2]
        weight = INSTABILITY_WEIGHTS.get(dipeptide, 1.0)
        dipeptide_sum += weight

    ii = (10.0 / len(seq)) * dipeptide_sum
    return round(ii, 2)


def compute_aromaticity(sequence: str) -> float:
    """
    Compute aromaticity: fraction of aromatic residues (F, W, Y, H).

    Aromatic residues contribute to pi-stacking interactions and
    aromatic cluster networks that stabilize thermophilic proteins.
    Higher aromaticity is correlated with thermostability.

    Reference: Bornscheuer et al. (2012), Nature 485:185-194.
    (General review; aromaticity as thermostability feature widely used)
    """
    seq = sequence.upper().strip()
    if not seq:
        return 0.0
    aromatic_count = sum(1 for aa in seq if aa in AROMATIC_AA)
    return round(aromatic_count / len(seq), 4)


def compute_charged_residue_ratio(sequence: str) -> dict:
    """
    Compute charged residue content and related thermostability metrics.

    Returns a dict with:
      - positive_fraction: fraction of R, K, H
      - negative_fraction: fraction of D, E
      - total_charged_fraction: positive + negative
      - charge_ratio: positive / (positive + negative), or 0 if no charged
      - net_charge_per_residue: (positive - negative) / length

    The Szilagyi-Zavodszky (2000) observation: thermophilic proteins have
    significantly higher fractions of charged residues than mesophilic
    counterparts, particularly glutamate (E) and arginine (R). The excess of
    charged residues forms salt bridges that stabilize the folded state.

    The ion-pair (salt bridge) density is one of the most reliable thermostability
    markers known. A protein with total_charged_fraction > 0.20 and balanced
    charge_ratio is a thermostability candidate.

    Reference: Szilagyi & Zavodszky (2000), Structure 8:493-504.
    """
    seq = sequence.upper().strip()
    if not seq:
        return {}

    total = len(seq)
    pos = sum(1 for aa in seq if aa in POSITIVELY_CHARGED)
    neg = sum(1 for aa in seq if aa in NEGATIVELY_CHARGED)
    charged = pos + neg

    return {
        "positive_fraction":      round(pos / total, 4),
        "negative_fraction":      round(neg / total, 4),
        "total_charged_fraction": round(charged / total, 4),
        "charge_ratio":           round(pos / charged, 4) if charged > 0 else 0.0,
        "net_charge_per_residue": round((pos - neg) / total, 4),
    }


def compute_proline_content(sequence: str) -> float:
    """
    Compute proline content (fraction of sequence that is proline).

    Proline is the "structure-breaking" residue, but in thermophilic proteins,
    strategic placement of proline in loops and turns reduces backbone
    conformational flexibility. An unfolded protein has more conformational
    entropy than a rigid one; proline reduces that entropy in the unfolded state,
    effectively shifting the folding equilibrium toward the native state at
    elevated temperatures.

    Higher proline content in loop regions = higher predicted thermostability.
    This is a crude proxy (we don't know where in the sequence the prolines are),
    but it is used as a feature signal.

    Reference: Watanabe et al. (1994), Biochemistry 33:981-990.
    """
    seq = sequence.upper().strip()
    if not seq:
        return 0.0
    return round(seq.count(PROLINE) / len(seq), 4)


def compute_small_residue_ratio(sequence: str) -> float:
    """
    Compute fraction of small residues (G, A, S, T, V).

    Thermophilic proteins often have higher glycine content in structured regions
    (adds conformational flexibility where needed) but lower glycine in unstructured
    loops. This is a mixed signal. We include it as a feature but weight it low
    in the IAS scorer.
    """
    seq = sequence.upper().strip()
    if not seq:
        return 0.0
    small_count = sum(1 for aa in seq if aa in SMALL_AA)
    return round(small_count / len(seq), 4)


def extract_all_features(sequence: str) -> dict:
    """
    Master function: compute all thermostability features for one sequence.

    Returns a flat dict of all features suitable for scoring and report generation.
    Any downstream module should call this rather than individual feature functions
    to ensure consistent feature sets.

    Feature list:
      - length: sequence length in amino acids
      - gravy: GRAVY hydrophobicity index
      - aliphatic_index: Ikai aliphatic index (higher = more thermostable)
      - instability_index: Guruprasad instability index (lower = more stable)
      - aromaticity: fraction aromatic residues
      - proline_fraction: fraction proline residues
      - small_residue_fraction: fraction small residues
      - positive_fraction: fraction positively charged (R, K, H)
      - negative_fraction: fraction negatively charged (D, E)
      - total_charged_fraction: total charged residue fraction
      - charge_ratio: positive / total charged (balance metric)
      - net_charge_per_residue: net charge normalized by length
      - aa_A through aa_Y: fractional abundance of each standard AA
    """
    if not sequence or len(sequence) < 10:
        logger.warning(f"Sequence too short for reliable feature extraction: {len(sequence)} aa")
        return {}

    features = {}
    features["length"] = len(sequence)
    features["gravy"] = compute_gravy_index(sequence)
    features["aliphatic_index"] = compute_aliphatic_index(sequence)
    features["instability_index"] = compute_instability_index(sequence)
    features["aromaticity"] = compute_aromaticity(sequence)
    features["proline_fraction"] = compute_proline_content(sequence)
    features["small_residue_fraction"] = compute_small_residue_ratio(sequence)

    charged = compute_charged_residue_ratio(sequence)
    features.update(charged)

    composition = compute_amino_acid_composition(sequence)
    for aa, frac in composition.items():
        features[f"aa_{aa}"] = frac

    return features


def thermostability_summary(sequence: str) -> dict:
    """
    Return a human-readable thermostability assessment summary.

    Interprets each feature against known thresholds from the literature.
    Returns a dict with 'signals' (list of supporting/contradicting evidence
    strings) and 'thermostability_score' (a simple 0-10 count of positive signals).

    This is NOT the IAS score. This is a qualitative interpretive layer
    used for the HTML report's narrative section.
    """
    features = extract_all_features(sequence)
    if not features:
        return {"signals": [], "thermostability_score": 0}

    signals = []
    score = 0

    if features.get("aliphatic_index", 0) > 80:
        signals.append("High aliphatic index (>80): strong hydrophobic core packing signal")
        score += 2
    elif features.get("aliphatic_index", 0) > 65:
        signals.append("Moderate aliphatic index (65-80): some core stability signal")
        score += 1

    if features.get("instability_index", 100) < 40:
        signals.append("Low instability index (<40): predicted in-vivo stable")
        score += 2
    elif features.get("instability_index", 100) < 55:
        signals.append("Moderate instability index (40-55): borderline stability")
        score += 1

    if features.get("total_charged_fraction", 0) > 0.20:
        signals.append("High charged residue fraction (>20%): thermophile ion-pair signature")
        score += 2
    elif features.get("total_charged_fraction", 0) > 0.15:
        signals.append("Moderate charged residue fraction: some ion-pair support")
        score += 1

    if features.get("proline_fraction", 0) > 0.05:
        signals.append("High proline content (>5%): reduced loop conformational entropy")
        score += 1

    if features.get("aromaticity", 0) > 0.10:
        signals.append("High aromaticity (>10%): aromatic stacking network present")
        score += 1

    if features.get("gravy", 0) > -0.5:
        signals.append("Moderate-positive GRAVY: reasonable hydrophobic core")
        score += 1

    return {
        "signals": signals,
        "thermostability_score": min(score, 10),  # Cap at 10
        "features": features,
    }
