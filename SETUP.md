# Archaeon Setup and Deployment Guide

Complete step-by-step instructions for getting Archaeon running, interpreting its output, and pushing to GitHub.

---

## Step 1: Get a Free NCBI API Key (Optional but Recommended)

Without an API key, NCBI rate-limits you to 3 requests per second. With one, you get 10 per second, which speeds up large runs significantly.

1. Go to [https://www.ncbi.nlm.nih.gov/account/](https://www.ncbi.nlm.nih.gov/account/)
2. Create a free NCBI account (email + password; no payment required)
3. After logging in, click your username in the top right, then "API Key Management"
4. Click "Create an API Key"
5. Copy the key that appears; you will paste it into your `.env` file in a moment

You also need an email address for NCBI's terms of service. It does not have to match your NCBI account email.

---

## Step 2: Set Up the Project Folder

Download or clone the Archaeon repository. If you are downloading the ZIP:

1. Download and extract the ZIP
2. You should have a folder called `archaeon` with this structure:

```
archaeon/
  src/
    data/
      ncbi.py
      mgnify.py
    analysis/
      features.py
      scorer.py
      structure.py
    visualization/
      report.py
  docs/
    lab_writeup.md
  archaeon.py
  requirements.txt
  .env.example
  .gitignore
  README.md
  SETUP.md
```

Open that folder in VS Code (File > Open Folder > select `archaeon`) or PyCharm (File > Open > select the `archaeon` folder).

---

## Step 3: Create a Virtual Environment

In the terminal inside the `archaeon` folder:

```bash
python -m venv venv
```

Then activate it:

```bash
# Windows
venv\Scripts\activate

# Mac/Linux
source venv/bin/activate
```

You should see `(venv)` appear in your terminal prompt. If you don't, the virtual environment is not active and the next step will install packages globally (which is fine but messier).

---

## Step 4: Install Dependencies

```bash
pip install -r requirements.txt
```

This installs biopython, requests, pandas, tqdm, colorama, and python-dotenv. It takes about 1-3 minutes. You will see a lot of output. Wait for it to finish and confirm no errors at the end.

If you see errors about a specific package, try:

```bash
pip install biopython requests pandas tqdm colorama python-dotenv
```

---

## Step 5: Create Your `.env` File

In the `archaeon` folder, copy `.env.example` to `.env`:

```bash
# Mac/Linux
cp .env.example .env

# Windows
copy .env.example .env
```

Open `.env` in any text editor. It should look like this:

```
NCBI_EMAIL=your_email@example.com
NCBI_API_KEY=
OUTPUT_DIR=data/outputs
```

Fill in your values:

```
NCBI_EMAIL=youremail@gmail.com
NCBI_API_KEY=abcdef1234567890abcdef1234567890abcd  (paste your key from Step 1, or leave blank)
OUTPUT_DIR=data/outputs
```

Save the file. Do not add spaces around the `=` sign. The `.env` file is in `.gitignore` and will never be uploaded to GitHub.

---

## Step 6: Run a Fast Test First

Before running the full pipeline, verify that the setup works with a minimal run:

```bash
python archaeon.py --max-per-biome 5 --top-n 10 --no-report
```

This fetches only 5 sequences per biome query (very fast) and skips the HTML report. You should see output like:

```
============================================================
  ARCHAEON: Extremophile Enzyme Discovery Pipeline
  The Xiu Lab | github.com/axshoe/archaeon
============================================================
  Run started: 2026-04-01 14:22:05
  Sources: ncbi | Max per biome: 5 | Top N: 10

============================================================
  Phase 1: Data Collection
============================================================
  Fetching from NCBI (max 5 per biome query)...
  Found 5 hits for 'thermophile[organism] AND enzyme[title]'
  ...
  Total unique candidates: 47

============================================================
  Phase 2: Sequence Feature Extraction
============================================================
  Features extracted for 47 candidates

============================================================
  Phase 3: IAS Scoring and Ranking
============================================================
  Top IAS: 73.4 (AAC69568.1)
  Mean IAS (top 10): 58.2
```

If you see errors about the NCBI API key, double-check that `.env` has no spaces around the `=` sign and that the file is named `.env` (not `.env.txt`).

---

## Step 7: Run the Standard Pipeline

```bash
python archaeon.py --sources ncbi --max-per-biome 20 --top-n 20
```

This takes about 2-4 minutes. Most of the time is NCBI API calls. At the end you should see:

```
============================================================
  PIPELINE COMPLETE
  Total time: 187.3s (3.1 minutes)
  Candidates evaluated: 156
  Top candidate: [accession] | IAS: 76.2 | Family: amylase
  Outputs: data/outputs/
============================================================
```

Open `data/outputs/archaeon_report.html` in your browser. This is your main output.

---

## Step 8: Add BLAST (Optional, Takes ~15-20 Minutes)

BLAST adds function annotation to the top candidates by comparing them to reference enzyme sequences on NCBI's servers. It makes the report significantly more informative but is slow.

```bash
python archaeon.py --max-per-biome 20 --top-n 20 --run-blast --blast-top-n 10
```

You can speed up repeat runs by caching the raw candidates from the first run and skipping data fetching:

```bash
# First run: save cache
python archaeon.py --max-per-biome 20 --save-cache data/candidates.json

# Second run: load cache, add BLAST
python archaeon.py --load-cache data/candidates.json --run-blast --blast-top-n 10
```

---

## Step 9: Add Structure Prediction (Optional, Takes ~5-10 Minutes)

Structure prediction runs ESMFold on the top candidates and embeds interactive 3D viewers in the HTML report.

```bash
python archaeon.py --max-per-biome 20 --run-structure --structure-top-n 5
```

The `--structure-top-n 5` flag means it only predicts structures for the top 5 candidates by IAS. Increase this number for more structures (each takes about 30-60 seconds).

---

## Step 10: Use Both Databases

To fetch from both NCBI and MGnify, which gives the broadest coverage of extremophile sequences:

```bash
python archaeon.py --sources both --max-per-biome 25 --top-n 20
```

Note: MGnify's protein API sometimes returns empty results if their servers are under load. If you see zero MGnify records, try again after a few minutes or run with `--sources ncbi` only.

---

## Step 11: Interpret the HTML Report

Open `data/outputs/archaeon_report.html` in any browser (Chrome or Firefox recommended for 3D visualization).

**Summary cards at the top:** Total candidates ranked, top IAS score, mean IAS score, number of structures predicted, top predicted enzyme family, top candidate organism.

**Score Distribution chart (left):** Histogram of all IAS scores in the candidate set. A good run will show a roughly normal distribution with a tail toward higher scores.

**Top 10 Candidates chart (center):** Horizontal bar chart of the top 10 candidates. Green bars are high-scoring (IAS > 70); amber bars are moderate; red bars are low.

**Aliphatic Index vs IAS scatter (right):** Confirms that aliphatic index is positively correlated with IAS (as expected given its weight in the formula). Candidates far above the trend line have strong BLAST matches in addition to structural features.

**Ranked candidate table:** Click any row to expand a detail panel showing the IAS component breakdown (thermostability / quality / BLAST) and individual sequence feature values.

**Structure viewers (if run-structure was enabled):** Interactive 3D molecular viewers. Colored by pLDDT confidence. Drag to rotate, scroll to zoom. High-confidence regions (dark blue) are the most reliable structural predictions.

**Methodology section at the bottom:** Documents the IAS formula, data sources, and limitations. Important to read before drawing conclusions.

---

## Step 12: Inspect the CSV

The CSV at `data/outputs/archaeon_candidates.csv` has one row per candidate with all scores and key features. Open in Excel or Google Sheets for custom analysis.

Key columns:
- `ias`: Final Industrial Applicability Score (0-100)
- `predicted_family`: Predicted enzyme family from BLAST
- `best_blast_identity`: % identity to best BLAST match
- `aliphatic_index`: Ikai aliphatic index (>80 = thermostable signal)
- `instability_index`: Guruprasad instability index (<40 = stable signal)
- `total_charged_fraction`: Fraction of charged residues (>0.20 = thermophile signal)
- `mean_plddt`: Mean ESMFold confidence score (only if structure prediction was run)

---

## Step 13: Push to GitHub

Create a new repository on GitHub called `archaeon`. Then in your terminal inside the `archaeon` folder:

```bash
git init
git add .
git commit -m "feat: Archaeon v1.0 - Extremophile Enzyme Discovery Pipeline"
git remote add origin https://github.com/YOUR_USERNAME/archaeon.git
git branch -M main
git push -u origin main
```

Your `.env` file will not be pushed (it's in `.gitignore`). Anyone cloning the repo will need to create their own `.env` from `.env.example`.

---

## Troubleshooting

**"NCBI rate limit exceeded" errors:** Add a free API key (Step 1) or increase the sleep time between requests by setting `--max-per-biome 10` (fewer requests per biome).

**"ESMFold API returned HTTP 503" errors:** EBI's ESMFold API is occasionally under maintenance. Try again after 10-15 minutes. If it persists, run without `--run-structure`.

**"No candidates retrieved":** Check that your `.env` file has the correct `NCBI_EMAIL` set. NCBI rejects anonymous requests. Also verify internet connectivity.

**BLAST takes too long:** BLAST on 10 candidates can take 15-30 minutes because each query hits NCBI's remote servers. Reduce `--blast-top-n` to 3-5 for faster runs. Or skip BLAST entirely; the pipeline works fine without it.

**Empty MGnify results:** MGnify's protein search API is in beta and occasionally returns empty results. Run with `--sources ncbi` if MGnify is giving you trouble. NCBI alone provides plenty of extremophile sequences.

**"ModuleNotFoundError: No module named 'src'":** Make sure you are running `python archaeon.py` from inside the `archaeon` folder, not from a parent directory.

---

## Minimum Viable Run (If Everything Else Fails)

```bash
python archaeon.py --sources ncbi --max-per-biome 5 --top-n 10
```

This fetches roughly 40-50 sequences, scores them, and generates a basic report in about 60 seconds. It won't have BLAST annotations or structure predictions, but it confirms the full pipeline runs end-to-end.

---

*Angie Xiu | The Xiu Lab | archaeon v1.0*
