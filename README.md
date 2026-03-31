# Archaeon

**Extremophile Enzyme Discovery Pipeline**

A Python CLI that mines NCBI and MGnify for protein sequences from extreme environments, scores them by predicted industrial applicability, optionally annotates them via BLAST, predicts 3D structures using ESMFold, and generates a self-contained interactive HTML report.

Built by [A. Xiu](https://github.com/axshoe) — [The Xiu Lab](https://thexiulab.org)

---

## Quickstart

```bash
git clone https://github.com/axshoe/archaeon
cd archaeon
pip install -r requirements.txt
cp .env.example .env   # add your NCBI email and API key
python archaeon.py --sources ncbi --max-per-biome 20 --top-n 20
python archaeon.py --sources ncbi --max-per-biome 20 --top-n 20
```

Open `data/outputs/archaeon_report.html` in any browser. No server needed.

---

## Easiest Access — No Setup Required

If you want to see Archaeon's output without installing anything, a pre-run report is hosted at:

**[axshoe.github.io/archaeon](https://axshoe.github.io/Archaeon/)**

This is a static HTML report generated from a real pipeline run (84 candidates, BLAST-annotated, 5 structures predicted via ESMFold). The 3D structure viewers are interactive — rotate, zoom, and pan with mouse/touch. No server, no login, no install.

---

## Installation

**Requirements:** Python 3.9+, internet connection, free NCBI account

```bash
# 1. Clone
git clone https://github.com/axshoe/archaeon
cd archaeon

# 2. Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate      # Mac/Linux
venv\Scripts\activate         # Windows PowerShell

# 3. Install dependencies
pip install -r requirements.txt

# 4. Configure
cp .env.example .env
# Edit .env — add NCBI_EMAIL and NCBI_API_KEY
# Get a free API key at: https://www.ncbi.nlm.nih.gov/account/
```

---

## Usage

### Fast run (no BLAST, no structure) — ~30 seconds
```bash
python archaeon.py --sources ncbi --max-per-biome 20 --top-n 20
```

### Standard run with structure prediction — ~3 minutes
```bash
python archaeon.py --sources ncbi --max-per-biome 20 --top-n 20 --run-structure --structure-top-n 5
```

### Full run with BLAST annotation — ~25–35 minutes
```bash
python archaeon.py --sources ncbi --max-per-biome 20 --top-n 20 --run-blast --blast-top-n 5 --run-structure --structure-top-n 5
```

### Both databases
```bash
python archaeon.py --sources both --max-per-biome 20 --top-n 20
```

---

## CLI Reference

| Flag | Default | Description |
|---|---|---|
| `--sources` | `ncbi` | Data sources: `ncbi`, `mgnify`, or `both` |
| `--max-per-biome` | `20` | Sequences retrieved per biome query |
| `--top-n` | `20` | Top N candidates to include in report |
| `--run-blast` | off | Enable BLAST annotation (adds 15–30 min) |
| `--blast-top-n` | `5` | Candidates to BLAST (top N by IAS) |
| `--run-structure` | off | Enable ESMFold structure prediction |
| `--structure-top-n` | `5` | Candidates to predict structures for |
| `--output-dir` | `data/outputs` | Output directory |
| `--no-report` | off | Skip HTML report generation |
| `--no-csv` | off | Skip CSV export |
| `--save-cache` | off | Cache retrieved sequences for reuse |
| `--load-cache` | off | Load sequences from cache |

---

## How It Works

**Phase 1 — Data retrieval.** Queries NCBI Entrez and/or MGnify for protein sequences from extreme environments: hydrothermal vents, hot springs, *Sulfolobus*, *Thermus thermophilus*, *Pyrococcus*, and others. Sequences are deduplicated across biome queries.

**Phase 2 — Feature extraction.** Computes five thermostability features per sequence from original literature formulas: GRAVY index (Kyte-Doolittle 1982), aliphatic index (Ikai 1980), instability index (Guruprasad 1990), charged residue fraction (Szilagyi & Zavodszky 2000), and proline content (Watanabe et al. 1996).

**Phase 3 — IAS scoring.** Ranks all candidates by the Industrial Applicability Score:

```
IAS = 0.40 × Thermostability + 0.20 × Quality + 0.40 × BLAST Identity
```

Thermostability combines the five features with weights derived from the literature. Quality is a piecewise linear function of sequence length (optimal 150–800 aa). BLAST Identity defaults to 50 (neutral) if BLAST is not run.

**Phase 4 — BLAST annotation (optional).** Submits top-N candidates to NCBI remote BLASTp against eight reference thermostable enzyme families: lipase, protease, amylase, cellulase, xylanase, laccase, peroxidase, isomerase. Updates IAS scores with real BLAST identity values.

**Phase 5 — Structure prediction (optional).** Submits top-N candidates to ESMFold's free public REST API. Parses per-residue pLDDT confidence scores from the B-factor column of the returned PDB file. Saves PDB files locally.

**Phase 6 — HTML report.** Generates a self-contained interactive report with IAS distribution charts, a ranked candidate table with collapsible feature detail, and 3D molecular viewers with pLDDT confidence coloring.

**Phase 7 — CSV export.** Exports all candidate data as a flat CSV for downstream analysis.

---

## Industrial Applicability Score (IAS)

The IAS is a composite ranking metric (0–100) designed for relative prioritization of enzyme candidates before experimental work. It is not a calibrated predictor of experimental thermostability — it is a structured summary of prior knowledge applied to sequence data.

| Component | Weight | Basis |
|---|---|---|
| Thermostability | 40% | 5 sequence features from extremophile literature |
| Sequence Quality | 20% | Length-based E. coli expression feasibility |
| BLAST Identity | 40% | Percent identity to reference thermostable enzymes |

Thermostability sub-weights:

| Feature | Weight | Paper |
|---|---|---|
| Aliphatic index | 30% | Ikai (1980) |
| Instability index (inverted) | 25% | Guruprasad et al. (1990) |
| Charged residue fraction | 25% | Szilagyi & Zavodszky (2000) |
| Proline content | 10% | Watanabe et al. (1996) |
| Aromaticity | 10% | — |

---

## Output Files

| File | Description |
|---|---|
| `data/outputs/archaeon_report.html` | Interactive HTML report — open in any browser |
| `data/outputs/archaeon_candidates.csv` | All candidates with features and scores |
| `data/outputs/structures/*.pdb` | PDB files for structure-predicted candidates |

---

## Limitations

- IAS weights are literature-informed priors, not fit to experimental thermostability data
- BLAST identity does not guarantee functional equivalence
- ESMFold pLDDT is a structural confidence score, not a thermostability predictor
- Sequence-based thermostability features are population-level generalizations; individual proteins frequently violate them
- All candidates require experimental validation before any industrial or commercial use

---

## References

| Paper | Used for |
|---|---|
| Altschul et al. (1990), *J Mol Biol* | BLAST algorithm |
| Bornscheuer et al. (2012), *Nature* | Industrial enzyme context |
| Guruprasad et al. (1990), *Protein Engineering* | Instability index |
| Handelsman (2004), *Microbiol Mol Biol Rev* | Metagenomics rationale |
| Ikai (1980), *J Biochem* | Aliphatic index |
| Kyte & Doolittle (1982), *J Mol Biol* | GRAVY / hydrophobicity scale |
| Lin et al. (2023), *Science* | ESMFold / pLDDT |
| Mitchell et al. (2020), *Nucleic Acids Res* | MGnify database |
| Sarmiento et al. (2015), *Front Bioeng Biotechnol* | Extremozyme applications |
| Sayers et al. (2022), *Nucleic Acids Res* | NCBI database |
| Szilagyi & Zavodszky (2000), *Structure* | Charged residue fraction |
| Watanabe et al. (1996), *Appl Environ Microbiol* | Proline content |

---

## Tech Stack

Python 3.9+ · BioPython · Requests · NCBI Entrez · MGnify REST API · ESMFold (Meta AI) · Chart.js · 3Dmol.js

---

*The Xiu Lab · 2026*
