# Archaeon

**Extremophile Enzyme Discovery Pipeline**

Archaeon is a Python CLI that mines public metagenomic databases (NCBI Entrez, MGnify) for novel enzyme candidates from extremophile organisms, predicts their 3D structure using ESMFold, scores them for industrial applicability using sequence-based thermostability features, and outputs a ranked candidate list with an interactive HTML report.

Everything runs from the command line. Everything is free. No GPU required.

---

## Why This Exists

Extremophile enzymes (from organisms living in volcanic vents, hot springs, hypersaline lakes) are enormously valuable for industrial biotechnology and drug manufacturing because they function under conditions that destroy normal enzymes. Finding them currently requires expensive wet-lab infrastructure or paid database subscriptions.

Archaeon compresses the computational discovery step into a fully reproducible open-source pipeline using only public databases (NCBI, MGnify) and free prediction tools (ESMFold via Meta AI's public REST API).

---

## Architecture

```
archaeon.py                     (CLI entry point, pipeline orchestration)
    |
    +-- src/data/
    |       ncbi.py             (NCBI Entrez API: esearch, efetch, BLASTp)
    |       mgnify.py           (MGnify JSON REST API: biome search, protein download)
    |
    +-- src/analysis/
    |       features.py         (Sequence feature extraction: GRAVY, AI, II, aromaticity, ...)
    |       scorer.py           (Industrial Applicability Score: weighted composite 0-100)
    |       structure.py        (ESMFold API: structure prediction, pLDDT parsing)
    |
    +-- src/visualization/
            report.py           (Self-contained HTML report: Chart.js + 3Dmol.js)

data/
    outputs/
        archaeon_report.html    (Main output: ranked candidates + visualizations)
        archaeon_candidates.csv (Flat CSV of all scored candidates)
        structures/             (PDB files from ESMFold predictions)
```

---

## Quick Start

```bash
# 1. Clone and enter the project
git clone https://github.com/axshoe/archaeon.git
cd archaeon

# 2. Create and activate a virtual environment
python -m venv venv
source venv/bin/activate       # Mac/Linux
venv\Scripts\activate          # Windows

# 3. Install dependencies
pip install -r requirements.txt

# 4. Configure credentials
cp .env.example .env
# Edit .env: add your NCBI email (required) and optional API key

# 5. Run the pipeline (fast mode: NCBI only, no BLAST, no structure prediction)
python archaeon.py

# 6. Open the report
open data/outputs/archaeon_report.html   # Mac
start data/outputs/archaeon_report.html  # Windows
```

---

## Usage

```
python archaeon.py [options]
```

| Flag | Default | Description |
|------|---------|-------------|
| `--sources` | `ncbi` | Data sources: `ncbi`, `mgnify`, or `both` |
| `--max-per-biome` | `20` | Max sequences per biome query |
| `--top-n` | `20` | Candidates to include in report |
| `--run-blast` | off | Run BLASTp on top candidates (slow: ~15 min) |
| `--blast-top-n` | `10` | How many top candidates to BLAST |
| `--run-structure` | off | Predict 3D structures via ESMFold |
| `--structure-top-n` | `5` | How many top candidates to predict structures for |
| `--load-cache` | None | Skip fetching; load candidates from JSON file |
| `--save-cache` | None | Save raw candidates to JSON after fetching |
| `--output-dir` | `data/outputs` | Output directory |
| `--verbose` | off | Debug logging |

### Example Runs

```bash
# Full pipeline: both databases, BLAST, structure prediction
python archaeon.py --sources both --max-per-biome 25 --run-blast --run-structure

# Fast exploration: NCBI only, skip BLAST
python archaeon.py --sources ncbi --max-per-biome 30 --top-n 25

# Load cached candidates (skip fetching, add BLAST)
python archaeon.py --load-cache data/candidates.json --run-blast --blast-top-n 5

# MGnify hydrothermal vent focus
python archaeon.py --sources mgnify --max-per-biome 50 --run-structure --structure-top-n 8
```

---

## Data Sources

| Source | Type | Access | Biomes Searched |
|--------|------|--------|-----------------|
| [NCBI Protein (nr)](https://www.ncbi.nlm.nih.gov/protein/) | Curated + environmental | Free API (email required) | Thermophile, hyperthermophile, hot spring, hydrothermal vent, halophile, acidophile, Sulfolobus, Thermus, Pyrococcus |
| [MGnify](https://www.ebi.ac.uk/metagenomics/) | Metagenomic exclusively | Free JSON REST API | Hydrothermal vents, hot springs, hypersaline, acid mine, deep subsurface, polar, volcanic |

**Why two sources:** NCBI's protein database skews toward cultured organisms with well-characterized genomes. MGnify is exclusively metagenomic; it captures sequences from organisms that have never been grown in a lab. The rarest extremophile diversity is uncultured. Using both maximizes candidate diversity.

---

## Industrial Applicability Score (IAS)

The IAS is a composite score (0-100) computed for each candidate:

```
IAS = 0.40 * Thermostability + 0.20 * Quality + 0.40 * BLAST
```

**Thermostability component** (40%):

| Feature | Weight | Rationale |
|---------|--------|-----------|
| Aliphatic index | 30% | Hydrophobic core packing; strongest single predictor |
| Instability index (inverted) | 25% | Dipeptide-based stability prediction |
| Total charged fraction | 25% | Ion pair / salt bridge density (thermophile signature) |
| Proline content | 10% | Loop rigidity; reduces unfolded state entropy |
| Aromaticity | 10% | Aromatic stacking network signal |

**Quality component** (20%): Sequence length penalty curve (optimal range: 150-800 aa for industrial expression in E. coli).

**BLAST component** (40%): Percent identity to the best-matching reference enzyme from 8 thermostable enzyme families (lipase, protease, amylase, cellulase, xylanase, laccase, peroxidase, isomerase). Higher identity = more confident function annotation.

---

## Thermostability Features: Formulas

**GRAVY Index** (Kyte & Doolittle, 1982):
```
GRAVY = (1/L) * sum_{i=1}^{L} H(aa_i)
```
Where H(aa) is the Kyte-Doolittle hydrophobicity scale. GRAVY > 0 = net hydrophobic.

**Aliphatic Index** (Ikai, 1980):
```
AI = X(A) + 2.9*X(V) + 3.9*(X(I) + X(L))
```
Where X(aa) = molar fraction of amino acid aa. AI > 80 typical for thermostable enzymes.

**Instability Index** (Guruprasad et al., 1990):
```
II = (10/L) * sum_{i=1}^{L-1} DIWV[aa_i][aa_{i+1}]
```
Where DIWV is the dipeptide instability weight table. II < 40 = predicted stable.

**Charged Residue Fraction** (Szilagyi & Zavodszky, 2000):
```
CRF = (|R| + |K| + |H| + |D| + |E|) / L
```
CRF > 0.20 is a thermophile signature; ion pairs stabilize the folded state.

---

## Structure Prediction

Structures are predicted using **ESMFold** (Lin et al., 2022, Science 379:1123-1130) via Meta AI's public REST API. Unlike AlphaFold2, ESMFold does not require a multiple sequence alignment, making it fast (seconds per sequence) and well-suited for metagenomic sequences from under-sampled biomes.

**pLDDT confidence scale:**

| Color | Range | Interpretation |
|-------|-------|----------------|
| Dark blue | ≥ 90 | Very high confidence |
| Light blue | 70-90 | Confident |
| Yellow | 50-70 | Low confidence (likely disordered) |
| Orange | < 50 | Very low confidence (do not interpret) |

Sequences >400 aa are truncated to the N-terminal 400 aa for API submission.

---

## Output Files

| File | Description |
|------|-------------|
| `data/outputs/archaeon_report.html` | Self-contained interactive HTML report (open in any browser) |
| `data/outputs/archaeon_candidates.csv` | Flat CSV with all scored candidates and features |
| `data/outputs/structures/*.pdb` | PDB structure files from ESMFold (one per predicted candidate) |
| `archaeon.log` | Run log with timing and API call details |

---

## Limitations

This is a computational first-pass tool. The IAS is a theoretical ranking, not an experimental thermostability measurement. Specific limitations:

- BLAST identity does not guarantee functional equivalence; homologs can have different substrate specificities
- IAS weights (40/20/40) are based on literature-informed priors, not calibrated against experimental data
- ESMFold pLDDT scores indicate structural confidence, not thermostability
- Sequences from environmental metagenomes may be incomplete ORFs (open reading frames), not full-length enzymes
- The aliphatic index threshold of >80 for thermostability is a population-level generalization; individual exceptions exist

All candidates identified by Archaeon require experimental validation before any industrial or commercial application.

---

## Project Context

Archaeon is the fourth project at [The Xiu Lab](https://github.com/axshoe), alongside:
- **DermEquity**: Medical AI fairness framework for skin cancer screening
- **NEXUS**: Quantitative market intelligence CLI (Node.js)
- **Stratum**: County-level Mobility Barrier Index with geospatial dashboard (Python)

---

## References

| Paper | Role in Pipeline |
|-------|-----------------|
| Kyte & Doolittle (1982). J. Mol. Biol. 157:105-132 | GRAVY hydrophobicity scale |
| Ikai (1980). J. Biochem. 88:1895-1898 | Aliphatic index formula |
| Guruprasad et al. (1990). Protein Engineering 4:155-161 | Instability index formula |
| Szilagyi & Zavodszky (2000). Structure 8:493-504 | Charged residue thermostability analysis |
| Lin et al. (2022). Science 379:1123-1130 | ESMFold structure prediction |
| Mitchell et al. (2020). Nucleic Acids Res. 48:D570-D578 | MGnify database |
| Sayers et al. (2022). Nucleic Acids Res. 50:D20-D26 | NCBI Entrez documentation |

---

## License

MIT. See LICENSE file.

---

*Angie Xiu | The Xiu Lab | Carmel, Indiana*
