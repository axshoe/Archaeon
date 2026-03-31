#!/usr/bin/env python3
"""
archaeon.py
-----------
Archaeon: Extremophile Enzyme Discovery Pipeline
Main CLI entry point.

Usage:
  python archaeon.py [options]

Examples:
  # Full pipeline: fetch, score, predict structures, generate report
  python archaeon.py --sources both --max-per-biome 20 --top-n 15 --run-blast --run-structure

  # Fast mode: NCBI only, no BLAST, no structure prediction
  python archaeon.py --sources ncbi --max-per-biome 30 --top-n 20

  # MGnify only with structure prediction on top 5
  python archaeon.py --sources mgnify --top-n 10 --run-structure --structure-top-n 5

  # Skip data fetching, load from saved CSV
  python archaeon.py --load-cache data/candidate_cache.json --run-blast --run-structure

Pipeline phases:
  Phase 1: Data collection (NCBI Entrez and/or MGnify API)
  Phase 2: Sequence feature extraction (thermostability features)
  Phase 3: IAS scoring and ranking
  Phase 4: BLAST similarity search on top candidates (optional, slow)
  Phase 5: ESMFold structure prediction on top candidates (optional)
  Phase 6: HTML report generation
  Phase 7: CSV export

Run time estimates (rough):
  Data collection (NCBI, 20/biome): ~2-3 minutes
  Data collection (MGnify, 25/biome): ~1-2 minutes
  Feature extraction (300 sequences): <1 second
  IAS scoring (300 sequences): <1 second
  BLAST (10 candidates): ~10-15 minutes (remote NCBI servers)
  ESMFold (10 candidates): ~2-3 minutes
  Total without BLAST: ~5 minutes
  Total with BLAST: ~20 minutes
"""

import os
import sys
import json
import csv
import logging
import argparse
from datetime import datetime
from pathlib import Path
from dotenv import load_dotenv

# Rich progress display
try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False
    # Fallback: simple print-based progress
    class tqdm:
        def __init__(self, iterable=None, **kwargs):
            self._iterable = iterable or []
            self._desc = kwargs.get("desc", "")
        def __enter__(self): return self
        def __exit__(self, *a): pass
        def __iter__(self): return iter(self._iterable)
        def set_description(self, d): print(f"  {d}")
        def update(self, n=1): pass

try:
    from colorama import Fore, Style, init as colorama_init
    colorama_init()
    COLORAMA_AVAILABLE = True
except ImportError:
    COLORAMA_AVAILABLE = False
    class Fore:
        CYAN = GREEN = YELLOW = RED = MAGENTA = WHITE = ""
    class Style:
        BRIGHT = RESET_ALL = ""

# Internal modules
from src.data.ncbi import configure_entrez, search_all_biomes, blast_against_enzyme_families
from src.data.mgnify import search_all_extreme_biomes
from src.analysis.features import extract_all_features
from src.analysis.scorer import rank_candidates
from src.analysis.structure import batch_predict_structures
from src.visualization.report import generate_report

# Load environment variables from .env file
load_dotenv()

# -------------------------------------------------------------------------
# Configure logging
# -------------------------------------------------------------------------
def setup_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler("archaeon.log", mode="a"),
        ]
    )


def cprint(msg: str, color: str = "", bold: bool = False) -> None:
    """Colorama-aware console print."""
    prefix = (Style.BRIGHT if bold else "") + color
    suffix = Style.RESET_ALL if (bold or color) else ""
    print(f"{prefix}{msg}{suffix}")


# -------------------------------------------------------------------------
# Phase banners for clean terminal output
# -------------------------------------------------------------------------
def phase_banner(phase: int, title: str) -> None:
    cprint(f"\n{'='*60}", Fore.CYAN)
    cprint(f"  Phase {phase}: {title}", Fore.CYAN, bold=True)
    cprint(f"{'='*60}", Fore.CYAN)


# -------------------------------------------------------------------------
# CLI argument parser
# -------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="archaeon",
        description="Archaeon: Extremophile Enzyme Discovery Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python archaeon.py                                # Default run (NCBI only, no BLAST)
  python archaeon.py --sources both --run-blast     # Full pipeline with BLAST
  python archaeon.py --run-structure --structure-top-n 5  # Structure prediction on top 5
  python archaeon.py --load-cache cache.json        # Skip fetching, use cached data
        """
    )

    # Data source options
    data_group = parser.add_argument_group("Data Sources")
    data_group.add_argument(
        "--sources", choices=["ncbi", "mgnify", "both"], default="ncbi",
        help="Which databases to fetch sequences from (default: ncbi)"
    )
    data_group.add_argument(
        "--max-per-biome", type=int, default=20,
        help="Max sequences per biome query (default: 20)"
    )
    data_group.add_argument(
        "--load-cache", type=str, default=None,
        help="Load candidates from a JSON cache file (skip fetching)"
    )
    data_group.add_argument(
        "--save-cache", type=str, default=None,
        help="Save raw candidates to a JSON file after fetching"
    )

    # Scoring options
    score_group = parser.add_argument_group("Scoring")
    score_group.add_argument(
        "--top-n", type=int, default=20,
        help="Number of top candidates to include in output (default: 20)"
    )
    score_group.add_argument(
        "--min-ias", type=float, default=0.0,
        help="Minimum IAS score threshold for inclusion in report (default: 0)"
    )

    # BLAST options
    blast_group = parser.add_argument_group("BLAST (optional, slow ~15min)")
    blast_group.add_argument(
        "--run-blast", action="store_true",
        help="Run BLAST similarity search on top candidates (adds ~15min)"
    )
    blast_group.add_argument(
        "--blast-top-n", type=int, default=10,
        help="Number of top candidates to run BLAST on (default: 10)"
    )

    # Structure prediction options
    struct_group = parser.add_argument_group("Structure Prediction (optional)")
    struct_group.add_argument(
        "--run-structure", action="store_true",
        help="Run ESMFold structure prediction on top candidates"
    )
    struct_group.add_argument(
        "--structure-top-n", type=int, default=5,
        help="Number of top candidates to predict structures for (default: 5)"
    )
    struct_group.add_argument(
        "--structure-output-dir", type=str, default="data/outputs/structures",
        help="Directory to save PDB files (default: data/outputs/structures)"
    )

    # Output options
    output_group = parser.add_argument_group("Output")
    output_group.add_argument(
        "--output-dir", type=str, default="data/outputs",
        help="Output directory for report and CSV (default: data/outputs)"
    )
    output_group.add_argument(
        "--report-name", type=str, default="archaeon_report.html",
        help="HTML report filename (default: archaeon_report.html)"
    )
    output_group.add_argument(
        "--no-report", action="store_true",
        help="Skip HTML report generation"
    )
    output_group.add_argument(
        "--no-csv", action="store_true",
        help="Skip CSV export"
    )

    # Misc
    parser.add_argument("--verbose", action="store_true", help="Verbose logging")
    parser.add_argument("--version", action="version", version="Archaeon 1.0.0")

    return parser


# -------------------------------------------------------------------------
# Phase 1: Data Collection
# -------------------------------------------------------------------------
def phase_fetch_data(args) -> list:
    """
    Fetch sequences from NCBI and/or MGnify.

    Returns a list of raw candidate dicts (before scoring).
    Each dict has at minimum: id, accession, description, sequence, length,
    organism, source_biome.
    """
    phase_banner(1, "Data Collection")

    # Check for NCBI credentials
    ncbi_email = os.getenv("NCBI_EMAIL")
    ncbi_api_key = os.getenv("NCBI_API_KEY")

    if not ncbi_email:
        cprint("WARNING: NCBI_EMAIL not set in .env. Using placeholder.", Fore.YELLOW)
        ncbi_email = "archaeon.pipeline@placeholder.edu"

    configure_entrez(ncbi_email, ncbi_api_key)

    all_candidates = []

    if args.sources in ("ncbi", "both"):
        cprint(f"  Fetching from NCBI (max {args.max_per_biome} per biome query)...", Fore.WHITE)
        ncbi_records = search_all_biomes(max_per_biome=args.max_per_biome)
        for rec in ncbi_records:
            rec["source"] = "NCBI"
        all_candidates.extend(ncbi_records)
        cprint(f"  {Fore.GREEN}NCBI: {len(ncbi_records)} sequences retrieved{Style.RESET_ALL}")

    if args.sources in ("mgnify", "both"):
        cprint(f"  Fetching from MGnify (max {args.max_per_biome} per biome)...", Fore.WHITE)
        mgnify_records = search_all_extreme_biomes(max_per_biome=args.max_per_biome)
        for rec in mgnify_records:
            if "source" not in rec:
                rec["source"] = "MGnify"
        all_candidates.extend(mgnify_records)
        cprint(f"  {Fore.GREEN}MGnify: {len(mgnify_records)} sequences retrieved{Style.RESET_ALL}")

    # Deduplicate by sequence (some sequences appear in both databases)
    seen_seqs = set()
    deduped = []
    for rec in all_candidates:
        seq_key = rec.get("sequence", "")[:50]  # Use first 50 aa as key
        if seq_key not in seen_seqs and seq_key:
            seen_seqs.add(seq_key)
            deduped.append(rec)

    cprint(f"\n  Total unique candidates: {len(deduped)} "
           f"(removed {len(all_candidates) - len(deduped)} duplicates)", Fore.CYAN, bold=True)
    return deduped


# -------------------------------------------------------------------------
# Phase 4: BLAST (optional)
# -------------------------------------------------------------------------
def phase_blast(candidates: list, top_n: int) -> dict:
    """
    Run BLAST on the top N candidates (pre-IAS-sorted).
    Returns dict mapping candidate id -> blast_hits list.
    """
    phase_banner(4, f"BLAST Similarity Search (top {top_n} candidates)")
    cprint("  Note: This uses NCBI remote BLAST and may take 15-30 minutes.", Fore.YELLOW)

    blast_results = {}
    targets = candidates[:top_n]

    for i, candidate in enumerate(targets):
        cid = candidate.get("id", f"c{i}")
        seq = candidate.get("sequence", "")
        cprint(f"  [{i+1}/{len(targets)}] BLASTing {cid} (length {len(seq)})...", Fore.WHITE)

        if not seq:
            continue

        hits = blast_against_enzyme_families(seq, top_n_families=3)
        blast_results[cid] = hits

        if hits:
            best = hits[0]
            cprint(f"    Best hit: {best['family']} ({best['identity_pct']}% identity)", Fore.GREEN)
        else:
            cprint(f"    No significant hits found", Fore.YELLOW)

    cprint(f"\n  BLAST complete for {len(blast_results)} candidates", Fore.CYAN, bold=True)
    return blast_results


# -------------------------------------------------------------------------
# Phase 5: Structure Prediction (optional)
# -------------------------------------------------------------------------
def phase_structure(candidates: list, args) -> dict:
    """
    Run ESMFold structure prediction on top candidates.
    Adds pdb_string, mean_plddt, quality_summary to each candidate.
    Returns dict mapping candidate_id -> structure result.
    """
    phase_banner(5, f"Structure Prediction (top {args.structure_top_n} candidates)")
    cprint("  Using ESMFold API (Meta AI Research) - free, no GPU needed.", Fore.WHITE)

    structure_results = batch_predict_structures(
        candidates=candidates,
        output_dir=args.structure_output_dir,
        max_candidates=args.structure_top_n,
    )

    successes = sum(1 for r in structure_results.values() if r.get("success"))
    cprint(f"\n  Structure prediction complete: {successes}/{len(structure_results)} succeeded",
           Fore.CYAN, bold=True)

    if successes > 0:
        for cid, res in structure_results.items():
            if res.get("success"):
                cprint(f"    {cid}: mean pLDDT = {res['mean_plddt']:.1f}", Fore.GREEN)

    return structure_results


# -------------------------------------------------------------------------
# CSV Export
# -------------------------------------------------------------------------
def export_csv(candidates: list, output_path: str) -> None:
    """
    Export ranked candidates to CSV.
    Flattens the 'features' nested dict into top-level columns.
    """
    if not candidates:
        return

    # Build flat rows
    flat_rows = []
    for cand in candidates:
        row = {
            "rank":                cand.get("rank", 0),
            "id":                  cand.get("id", ""),
            "accession":           cand.get("accession", ""),
            "description":         cand.get("description", ""),
            "organism":            cand.get("organism", ""),
            "source_biome":        cand.get("source_biome", ""),
            "source":              cand.get("source", ""),
            "length":              cand.get("length", 0),
            "ias":                 cand.get("ias", 0),
            "thermostability_score": cand.get("thermostability_score", 0),
            "quality_score":       cand.get("quality_score", 0),
            "blast_score":         cand.get("blast_score", 0),
            "predicted_family":    cand.get("predicted_family", ""),
            "best_blast_identity": cand.get("best_blast_identity", 0),
            "mean_plddt":          cand.get("mean_plddt", 0),
        }
        # Add key sequence features
        features = cand.get("features", {})
        for key in ["aliphatic_index", "instability_index", "gravy", "aromaticity",
                    "proline_fraction", "total_charged_fraction", "charge_ratio"]:
            row[key] = features.get(key, "")
        flat_rows.append(row)

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(flat_rows[0].keys()))
        writer.writeheader()
        writer.writerows(flat_rows)

    cprint(f"  CSV exported to {output_path}", Fore.GREEN)


# -------------------------------------------------------------------------
# Main
# -------------------------------------------------------------------------
def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    setup_logging(args.verbose)

    start_time = datetime.now()

    cprint("\n" + "="*60, Fore.CYAN)
    cprint("  ARCHAEON: Extremophile Enzyme Discovery Pipeline", Fore.CYAN, bold=True)
    cprint("  The Xiu Lab | github.com/axshoe/archaeon", Fore.CYAN)
    cprint("="*60, Fore.CYAN)
    cprint(f"  Run started: {start_time.strftime('%Y-%m-%d %H:%M:%S')}", Fore.WHITE)
    cprint(f"  Sources: {args.sources} | Max per biome: {args.max_per_biome} | "
           f"Top N: {args.top_n}", Fore.WHITE)

    # -----------------------------------------------------------------------
    # Phase 1: Data Collection (or load from cache)
    # -----------------------------------------------------------------------
    if args.load_cache:
        phase_banner(1, "Loading Cached Data")
        with open(args.load_cache, "r") as f:
            all_candidates = json.load(f)
        cprint(f"  Loaded {len(all_candidates)} candidates from {args.load_cache}", Fore.GREEN)
    else:
        all_candidates = phase_fetch_data(args)

        if args.save_cache:
            os.makedirs(os.path.dirname(args.save_cache) or ".", exist_ok=True)
            with open(args.save_cache, "w") as f:
                json.dump(all_candidates, f, indent=2)
            cprint(f"  Cached {len(all_candidates)} candidates to {args.save_cache}", Fore.GREEN)

    if not all_candidates:
        cprint("\nERROR: No candidates retrieved. Check API credentials and network.", Fore.RED, bold=True)
        return 1

    # -----------------------------------------------------------------------
    # Phase 2: Feature Extraction
    # -----------------------------------------------------------------------
    phase_banner(2, "Sequence Feature Extraction")
    cprint(f"  Extracting thermostability features for {len(all_candidates)} sequences...", Fore.WHITE)

    for cand in all_candidates:
        seq = cand.get("sequence", "")
        if seq and "features" not in cand:
            cand["features"] = extract_all_features(seq)

    cprint(f"  {Fore.GREEN}Features extracted for {len(all_candidates)} candidates{Style.RESET_ALL}")

    # -----------------------------------------------------------------------
    # Phase 3: Scoring and Ranking
    # -----------------------------------------------------------------------
    phase_banner(3, "IAS Scoring and Ranking")
    cprint(f"  Computing Industrial Applicability Scores...", Fore.WHITE)

    ranked = rank_candidates(
        records=all_candidates,
        blast_results=None,  # BLAST not run yet; will re-rank after if BLAST enabled
        top_n=args.top_n * 3,  # Get extra candidates for BLAST selection
    )

    # Add rank numbers
    for i, cand in enumerate(ranked):
        cand["rank"] = i + 1

    cprint(f"  Top IAS: {ranked[0]['ias']:.1f} ({ranked[0].get('id','?')})", Fore.GREEN)
    cprint(f"  Mean IAS (top {len(ranked)}): {sum(c['ias'] for c in ranked)/len(ranked):.1f}", Fore.WHITE)

    # -----------------------------------------------------------------------
    # Phase 4: BLAST (optional)
    # -----------------------------------------------------------------------
    blast_results = {}
    if args.run_blast:
        blast_results = phase_blast(ranked, top_n=args.blast_top_n)

        # Re-rank with BLAST data incorporated
        cprint("\n  Re-ranking with BLAST scores...", Fore.WHITE)
        ranked = rank_candidates(
            records=all_candidates,
            blast_results=blast_results,
            top_n=args.top_n,
        )
        for i, cand in enumerate(ranked):
            cand["rank"] = i + 1
    else:
        cprint("\nPhase 4: BLAST - SKIPPED (use --run-blast to enable)", Fore.YELLOW)
        ranked = ranked[:args.top_n]

    # -----------------------------------------------------------------------
    # Phase 5: Structure Prediction (optional)
    # -----------------------------------------------------------------------
    structure_results = {}
    if args.run_structure:
        structure_results = phase_structure(ranked, args)
    else:
        cprint("\nPhase 5: Structure Prediction - SKIPPED (use --run-structure to enable)", Fore.YELLOW)

    # -----------------------------------------------------------------------
    # Phase 6: HTML Report
    # -----------------------------------------------------------------------
    if not args.no_report:
        phase_banner(6, "HTML Report Generation")
        report_path = os.path.join(args.output_dir, args.report_name)
        run_metadata = {
            "timestamp":          start_time.strftime("%Y-%m-%d %H:%M:%S"),
            "n_total_candidates": len(all_candidates),
            "sources":            args.sources.upper(),
            "version":            "1.0.0",
        }

        # Merge structure results into ranked candidates
        if structure_results:
            for cand in ranked:
                cid = cand.get("id", "")
                if cid in structure_results and structure_results[cid].get("success"):
                    cand["pdb_string"]     = structure_results[cid].get("pdb_string", "")
                    cand["mean_plddt"]     = structure_results[cid].get("mean_plddt", 0)
                    cand["quality_summary"] = structure_results[cid].get("quality_summary", {})

        generate_report(
            candidates=ranked,
            structure_results=structure_results,
            output_path=report_path,
            run_metadata=run_metadata,
        )
        cprint(f"  {Fore.GREEN}Report written: {report_path}{Style.RESET_ALL}")

    # -----------------------------------------------------------------------
    # Phase 7: CSV Export
    # -----------------------------------------------------------------------
    if not args.no_csv:
        phase_banner(7, "CSV Export")
        csv_path = os.path.join(args.output_dir, "archaeon_candidates.csv")
        export_csv(ranked, csv_path)

    # -----------------------------------------------------------------------
    # Summary
    # -----------------------------------------------------------------------
    elapsed = (datetime.now() - start_time).total_seconds()
    cprint("\n" + "="*60, Fore.GREEN)
    cprint("  PIPELINE COMPLETE", Fore.GREEN, bold=True)
    cprint(f"  Total time: {elapsed:.1f}s ({elapsed/60:.1f} minutes)", Fore.GREEN)
    cprint(f"  Candidates evaluated: {len(all_candidates)}", Fore.GREEN)
    cprint(f"  Top candidate: {ranked[0].get('id','?')} | IAS: {ranked[0]['ias']:.1f} | "
           f"Family: {ranked[0].get('predicted_family','unknown')}", Fore.GREEN)
    cprint(f"  Outputs: {args.output_dir}/", Fore.GREEN)
    cprint("="*60, Fore.GREEN)

    return 0


if __name__ == "__main__":
    sys.exit(main())
