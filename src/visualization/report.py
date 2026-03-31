"""
src/visualization/report.py
----------------------------
HTML report generator for Archaeon pipeline output.

Generates a self-contained HTML file with:
  1. Ranked candidate table (sortable, with IAS score and component breakdown)
  2. Sequence feature radar/bar charts per candidate (Chart.js, CDN-loaded)
  3. 3D structure viewer for each candidate with a predicted structure (3Dmol.js)
  4. IAS component distribution summary
  5. Data provenance and methodology notes

Design philosophy:
  Self-contained output. The HTML file loads Chart.js and 3Dmol.js from CDN,
  embeds all data inline as JavaScript variables, and requires no server to
  view. You can email this file, commit it to GitHub Pages, or open it in
  any browser without dependencies.

  The 3D structure viewer uses py3Dmol's JavaScript output (3Dmol.js library).
  PDB data is embedded directly in the HTML as a JavaScript string and loaded
  with viewer.addModel(pdbString, 'pdb'). This avoids file path issues.

  The report is intentionally over-documented. Future readers (including me
  reading this in six months wondering what I did) should understand what
  every number means without reading the source code.
"""

import os
import json
import math
import logging
from datetime import datetime
from typing import Optional
from pathlib import Path

logger = logging.getLogger(__name__)

# Color scheme: deep teal science lab aesthetic
# Inspired by genome browser color conventions + industrial design
COLORS = {
    "primary":    "#1a6b72",  # Deep teal
    "secondary":  "#2c9e8e",  # Medium teal
    "accent":     "#e8a838",  # Amber accent
    "high_score": "#2ecc71",  # Green for high IAS
    "mid_score":  "#f39c12",  # Amber for mid IAS
    "low_score":  "#e74c3c",  # Red for low IAS
    "bg":         "#0d1117",  # Dark background (GitHub dark mode inspired)
    "bg_card":    "#161b22",  # Card background
    "text":       "#c9d1d9",  # Primary text
    "text_muted": "#8b949e",  # Secondary text
    "border":     "#30363d",  # Border color
}


def _ias_color(ias: float) -> str:
    """Return a hex color based on IAS score for visual ranking."""
    if ias >= 70:
        return COLORS["high_score"]
    elif ias >= 45:
        return COLORS["mid_score"]
    else:
        return COLORS["low_score"]


def _format_number(val, decimals: int = 2) -> str:
    """Format a number for display; return 'N/A' for None/empty."""
    if val is None:
        return "N/A"
    try:
        return f"{float(val):.{decimals}f}"
    except (TypeError, ValueError):
        return str(val)


def _build_candidate_row(candidate: dict, rank: int) -> str:
    """
    Build a single HTML table row for a candidate.
    The row is collapsible: clicking reveals the feature breakdown.
    """
    ias = candidate.get("ias", 0)
    cid = candidate.get("id", f"c{rank}")
    desc = candidate.get("description", "")[:60]
    organism = candidate.get("organism", "unknown")[:40]
    biome = candidate.get("source_biome", "unknown")[:40]
    family = candidate.get("predicted_family", "unknown")
    identity = candidate.get("best_blast_identity", 0)
    length = candidate.get("length", 0)
    ai = candidate.get("features", {}).get("aliphatic_index", 0)
    ii = candidate.get("features", {}).get("instability_index", 0)
    tcf = candidate.get("features", {}).get("total_charged_fraction", 0)

    thermo_score = candidate.get("thermostability_score", 0)
    quality_score = candidate.get("quality_score", 0)
    blast_score = candidate.get("blast_score", 0)

    color = _ias_color(ias)

    return f"""
    <tr class="candidate-row" onclick="toggleDetail('detail-{rank}')">
      <td class="rank-cell">#{rank}</td>
      <td><code class="accession">{cid}</code></td>
      <td class="desc-cell">{desc}</td>
      <td>{organism}</td>
      <td><span class="biome-tag">{biome}</span></td>
      <td><span class="family-tag">{family}</span></td>
      <td class="score-cell" style="color:{color};font-weight:bold;">{_format_number(ias)}</td>
<td>{f"{identity:.2f}%" if identity > 0 else "&mdash;"}</td>      <td>{length}</td>
    </tr>
    <tr class="detail-row" id="detail-{rank}" style="display:none;">
      <td colspan="9">
        <div class="detail-panel">
          <div class="detail-grid">
            <div class="detail-section">
              <h4>IAS Component Breakdown</h4>
              <div class="component-bar">
                <label>Thermostability (40%)</label>
                <div class="bar-track">
                  <div class="bar-fill thermo-bar" style="width:{min(thermo_score,100)}%"></div>
                </div>
                <span>{_format_number(thermo_score)}</span>
              </div>
              <div class="component-bar">
                <label>Quality (20%)</label>
                <div class="bar-track">
                  <div class="bar-fill quality-bar" style="width:{min(quality_score,100)}%"></div>
                </div>
                <span>{_format_number(quality_score)}</span>
              </div>
              <div class="component-bar">
                <label>BLAST Identity (40%)</label>
                <div class="bar-track">
                  <div class="bar-fill blast-bar" style="width:{min(blast_score,100)}%"></div>
                </div>
                <span>{_format_number(blast_score)}</span>
              </div>
            </div>
            <div class="detail-section">
              <h4>Sequence Features</h4>
              <table class="feature-table">
                <tr><td>Aliphatic Index</td><td>{_format_number(ai)}</td></tr>
                <tr><td>Instability Index</td><td>{_format_number(ii)}</td></tr>
                <tr><td>Total Charged Fraction</td><td>{_format_number(tcf, 3)}</td></tr>
                <tr><td>GRAVY</td><td>{_format_number(candidate.get("features",{}).get("gravy",0))}</td></tr>
                <tr><td>Aromaticity</td><td>{_format_number(candidate.get("features",{}).get("aromaticity",0),3)}</td></tr>
                <tr><td>Proline Fraction</td><td>{_format_number(candidate.get("features",{}).get("proline_fraction",0),3)}</td></tr>
              </table>
            </div>
          </div>
          {"<div class='signals-section'><h4>Thermostability Signals</h4><ul>" + "".join(f"<li>{s}</li>" for s in candidate.get("thermostability_signals",[])) + "</ul></div>" if candidate.get("thermostability_signals") else ""}
        </div>
      </td>
    </tr>
    """


def _build_structure_viewer(candidate: dict, rank: int) -> str:
    pdb_string = candidate.get("pdb_string", "")
    if not pdb_string:
        return ""

    mean_plddt = candidate.get("mean_plddt", 0)
    qs = candidate.get("quality_summary", {})
    cid = candidate.get("id", f"c{rank}")

    # Rescale B-factors from 0-1 to 0-100 if ESMFold returned fractional pLDDT
    def _rescale_bfactors(pdb):
        lines = pdb.split('\n')
        out = []
        needs_rescale = None
        for line in lines:
            if line.startswith('ATOM') and needs_rescale is None:
                try:
                    b = float(line[60:66])
                    needs_rescale = b <= 1.0
                except (ValueError, IndexError):
                    pass
            if line.startswith('ATOM') and needs_rescale:
                try:
                    b = float(line[60:66])
                    scaled = f"{b * 100:6.2f}"
                    line = line[:60] + scaled + line[66:]
                except (ValueError, IndexError):
                    pass
            out.append(line)
        return '\n'.join(out)

    pdb_string = _rescale_bfactors(pdb_string)
    pdb_escaped = pdb_string.replace("\\", "\\\\").replace("`", "\\`").replace("${", "\\${")

    return f"""
    <div class="structure-card">
      <h3>Structure: {cid} (Rank #{rank})</h3>
      <div class="structure-meta">
        Mean pLDDT: <strong>{mean_plddt:.1f}</strong> |
        Very High: {qs.get('very_high',0)*100:.0f}% |
        High: {qs.get('high',0)*100:.0f}% |
        Low: {qs.get('low',0)*100:.0f}% |
        Very Low: {qs.get('very_low',0)*100:.0f}%
      </div>
      <div id="viewer-{rank}" class="structure-viewer"></div>
      <div class="plddt-legend">
        <span class="plddt-swatch" style="background:#0053D6;"></span> Very High (&ge;90)
        <span class="plddt-swatch" style="background:#65CBF3;"></span> High (70-90)
        <span class="plddt-swatch" style="background:#FFDB13;"></span> Low (50-70)
        <span class="plddt-swatch" style="background:#FF7D45;"></span> Very Low (&lt;50)
      </div>
      <script type="text/pdb-data" id="pdb-data-{rank}">{pdb_escaped}</script>
    </div>
    """


def _build_summary_charts_js(candidates: list) -> str:
    """
    Generate Chart.js JavaScript for summary visualizations.

    Creates:
      1. IAS score distribution histogram (all candidates)
      2. Top-10 candidates horizontal bar chart
      3. Feature correlation scatter plot (aliphatic index vs IAS)
    """
    # Extract data for charts
    ias_scores = [c.get("ias", 0) for c in candidates]
    ids = [c.get("id", f"c{i}")[:15] for i, c in enumerate(candidates[:10])]
    top10_ias = [c.get("ias", 0) for c in candidates[:10]]
    ai_values = [c.get("features", {}).get("aliphatic_index", 0) for c in candidates[:20]]
    ias_20 = [c.get("ias", 0) for c in candidates[:20]]
    scatter_data = [{"x": float(ai), "y": float(ias)} for ai, ias in zip(ai_values, ias_20)]

    return f"""
    // Chart 1: IAS Distribution Histogram
    (function() {{
        var scores = {json.dumps(ias_scores)};
        var bins = Array(10).fill(0);
        scores.forEach(function(s) {{
            var bin = Math.min(Math.floor(s / 10), 9);
            bins[bin]++;
        }});
        var ctx = document.getElementById('chart-distribution').getContext('2d');
        new Chart(ctx, {{
            type: 'bar',
            data: {{
                labels: ['0-10','10-20','20-30','30-40','40-50','50-60','60-70','70-80','80-90','90-100'],
                datasets: [{{
                    label: 'Number of Candidates',
                    data: bins,
                    backgroundColor: 'rgba(44, 158, 142, 0.7)',
                    borderColor: '#2c9e8e',
                    borderWidth: 1,
                }}]
            }},
            options: {{
                responsive: true,
                maintainAspectRatio: false,
                plugins: {{ legend: {{ display: false }}, title: {{ display: true, text: 'IAS Score Distribution', color: '#c9d1d9' }} }},
                scales: {{
                    x: {{ ticks: {{ color: '#8b949e' }}, grid: {{ color: '#30363d' }} }},
                    y: {{ ticks: {{ color: '#8b949e' }}, grid: {{ color: '#30363d' }} }}
                }}
            }}
        }});
    }})();

    // Chart 2: Top-10 Candidates
    (function() {{
        var ctx = document.getElementById('chart-top10').getContext('2d');
        var colors = {json.dumps([_ias_color(s) for s in top10_ias])};
        new Chart(ctx, {{
            type: 'bar',
            data: {{
                labels: {json.dumps(ids)},
                datasets: [{{
                    label: 'IAS Score',
                    data: {json.dumps(top10_ias)},
                    backgroundColor: colors,
                    borderColor: colors,
                    borderWidth: 1,
                }}]
            }},
            options: {{
                indexAxis: 'y',
                responsive: true,
                maintainAspectRatio: false,
                plugins: {{ legend: {{ display: false }}, title: {{ display: true, text: 'Top 10 Candidates by IAS', color: '#c9d1d9' }} }},
                scales: {{
                    x: {{ min: 0, max: 100, ticks: {{ color: '#8b949e' }}, grid: {{ color: '#30363d' }} }},
                    y: {{ ticks: {{ color: '#8b949e', font: {{ size: 11 }} }}, grid: {{ color: '#30363d' }} }}
                }}
            }}
        }});
    }})();

    // Chart 3: Aliphatic Index vs IAS scatter
    (function() {{
        var ctx = document.getElementById('chart-scatter').getContext('2d');
        var data = {json.dumps(scatter_data)};
        new Chart(ctx, {{
            type: 'scatter',
            data: {{
                datasets: [{{
                    label: 'Aliphatic Index vs IAS',
                    data: data,
                    backgroundColor: 'rgba(232, 168, 56, 0.6)',
                    borderColor: '#e8a838',
                    pointRadius: 5,
                }}]
            }},
            options: {{
                responsive: true,
                maintainAspectRatio: false,
                plugins: {{ legend: {{ display: false }}, title: {{ display: true, text: 'Aliphatic Index vs IAS (Top 20)', color: '#c9d1d9' }} }},
                scales: {{
                    x: {{ title: {{ display: true, text: 'Aliphatic Index', color: '#8b949e' }}, ticks: {{ color: '#8b949e' }}, grid: {{ color: '#30363d' }} }},
                    y: {{ title: {{ display: true, text: 'IAS Score', color: '#8b949e' }}, ticks: {{ color: '#8b949e' }}, grid: {{ color: '#30363d' }} }}
                }}
            }}
        }});
    }})();
    """


def generate_report(
    candidates: list,
    structure_results: Optional[dict] = None,
    output_path: str = "archaeon_report.html",
    run_metadata: Optional[dict] = None,
) -> str:
    """
    Generate the full HTML report.

    Args:
        candidates:       Ranked list from scorer.rank_candidates()
        structure_results: Dict from structure.batch_predict_structures()
                           mapping candidate_id -> structure result.
                           None if structure prediction was skipped.
        output_path:      File path for the output HTML file
        run_metadata:     Dict with run info: {timestamp, n_sources, biomes_searched, etc.}

    Returns:
        Path to the generated HTML file.
    """
    if not candidates:
        logger.error("No candidates to report.")
        return ""

    # Merge structure results into candidate dicts if available
    if structure_results:
        for cand in candidates:
            cid = cand.get("id", "")
            if cid in structure_results and structure_results[cid].get("success"):
                cand["pdb_string"]    = structure_results[cid].get("pdb_string", "")
                cand["mean_plddt"]    = structure_results[cid].get("mean_plddt", 0)
                cand["quality_summary"] = structure_results[cid].get("quality_summary", {})

    meta = run_metadata or {}
    timestamp = meta.get("timestamp", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    n_total = meta.get("n_total_candidates", len(candidates))
    sources = meta.get("sources", "NCBI Protein, MGnify")
    version = meta.get("version", "1.0.0")

    # Build candidate table rows
    table_rows = ""
    for i, cand in enumerate(candidates):
        table_rows += _build_candidate_row(cand, i + 1)

    # Build structure viewers (only for candidates with structure data)
    structure_viewers = ""
    structure_count = 0
    for i, cand in enumerate(candidates):
        if cand.get("pdb_string"):
            structure_viewers += _build_structure_viewer(cand, i + 1)
            structure_count += 1

    # Summary statistics
    mean_ias = sum(c.get("ias", 0) for c in candidates) / max(len(candidates), 1)
    max_ias = max(c.get("ias", 0) for c in candidates) if candidates else 0
    top_family = candidates[0].get("predicted_family", "unknown") if candidates else "N/A"
    top_organism = candidates[0].get("organism", "unknown") if candidates else "N/A"

    charts_js = _build_summary_charts_js(candidates)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Archaeon | Extremophile Enzyme Discovery Report</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.3/3Dmol-min.js"></script>
<style>
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ font-family: 'Segoe UI', system-ui, sans-serif; background: {COLORS["bg"]}; color: {COLORS["text"]}; line-height: 1.6; }}
  .header {{ background: linear-gradient(135deg, {COLORS["primary"]}, #0d2c2f); padding: 40px 60px; border-bottom: 1px solid {COLORS["border"]}; }}
  .header h1 {{ font-size: 2.2em; font-weight: 700; letter-spacing: -0.5px; }}
  .header h1 span {{ color: {COLORS["secondary"]}; }}
  .header p {{ color: {COLORS["text_muted"]}; margin-top: 8px; font-size: 0.95em; }}
  .container {{ max-width: 1400px; margin: 0 auto; padding: 30px 40px; }}
  .stat-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 16px; margin: 24px 0; }}
  .stat-card {{ background: {COLORS["bg_card"]}; border: 1px solid {COLORS["border"]}; border-radius: 8px; padding: 20px; text-align: center; }}
  .stat-card .stat-value {{ font-size: 2em; font-weight: 700; color: {COLORS["secondary"]}; }}
  .stat-card .stat-label {{ color: {COLORS["text_muted"]}; font-size: 0.85em; margin-top: 4px; }}
  .section-title {{ font-size: 1.4em; font-weight: 600; color: {COLORS["secondary"]}; margin: 32px 0 16px; padding-bottom: 8px; border-bottom: 1px solid {COLORS["border"]}; }}
  .charts-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(380px, 1fr)); gap: 24px; margin: 20px 0; }}
  .chart-card {{ background: {COLORS["bg_card"]}; border: 1px solid {COLORS["border"]}; border-radius: 8px; padding: 20px; }}
  .chart-card canvas {{ max-height: 280px; }}
  table.candidates-table {{ width: 100%; border-collapse: collapse; font-size: 0.9em; }}
  table.candidates-table th {{ background: {COLORS["primary"]}; color: white; padding: 10px 12px; text-align: left; font-weight: 600; position: sticky; top: 0; }}
  table.candidates-table td {{ padding: 10px 12px; border-bottom: 1px solid {COLORS["border"]}; }}
  tr.candidate-row:hover {{ background: rgba(44,158,142,0.08); cursor: pointer; }}
  .rank-cell {{ font-weight: 700; color: {COLORS["text_muted"]}; }}
  code.accession {{ font-family: monospace; font-size: 0.9em; background: rgba(255,255,255,0.05); padding: 2px 6px; border-radius: 3px; }}
  .biome-tag, .family-tag {{ background: rgba(26,107,114,0.3); border: 1px solid {COLORS["primary"]}; border-radius: 12px; padding: 2px 8px; font-size: 0.82em; white-space: nowrap; }}
  .family-tag {{ background: rgba(232,168,56,0.2); border-color: {COLORS["accent"]}; }}
  .score-cell {{ font-size: 1.1em; }}
  .detail-panel {{ background: rgba(22,27,34,0.9); border: 1px solid {COLORS["border"]}; border-radius: 6px; padding: 20px; margin: 4px 0; }}
  .detail-grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 24px; }}
  .detail-section h4 {{ color: {COLORS["secondary"]}; margin-bottom: 12px; font-size: 0.95em; text-transform: uppercase; letter-spacing: 0.5px; }}
  .component-bar {{ display: flex; align-items: center; gap: 10px; margin: 8px 0; font-size: 0.88em; }}
  .component-bar label {{ width: 160px; color: {COLORS["text_muted"]}; flex-shrink: 0; }}
  .bar-track {{ flex: 1; height: 8px; background: {COLORS["border"]}; border-radius: 4px; overflow: hidden; }}
  .bar-fill {{ height: 100%; border-radius: 4px; transition: width 0.3s; }}
  .thermo-bar {{ background: {COLORS["secondary"]}; }}
  .quality-bar {{ background: #9b59b6; }}
  .blast-bar {{ background: {COLORS["accent"]}; }}
  .component-bar span {{ width: 40px; text-align: right; font-weight: 600; font-size: 0.9em; }}
  .feature-table {{ width: 100%; font-size: 0.88em; }}
  .feature-table td {{ padding: 4px 8px; border-bottom: 1px solid rgba(48,54,61,0.5); }}
  .feature-table td:last-child {{ text-align: right; font-weight: 600; color: {COLORS["secondary"]}; }}
  .signals-section {{ margin-top: 16px; }}
  .signals-section ul {{ list-style: none; font-size: 0.85em; color: {COLORS["text_muted"]}; }}
  .signals-section li {{ padding: 3px 0; }}
  .signals-section li::before {{ content: "✓ "; color: {COLORS["secondary"]}; }}
  .structure-card {{ background: {COLORS["bg_card"]}; border: 1px solid {COLORS["border"]}; border-radius: 8px; padding: 20px; margin: 16px 0; }}
  .structure-card h3 {{ color: {COLORS["secondary"]}; margin-bottom: 8px; }}
  .structure-meta {{ font-size: 0.88em; color: {COLORS["text_muted"]}; margin-bottom: 12px; }}
  .structure-viewer {{ width: 100%; height: 400px; border-radius: 6px; border: 1px solid {COLORS["border"]}; position: relative; overflow: hidden; }}  .plddt-legend {{ display: flex; gap: 20px; margin-top: 10px; font-size: 0.82em; color: {COLORS["text_muted"]}; align-items: center; }}
  .plddt-swatch {{ display: inline-block; width: 12px; height: 12px; border-radius: 2px; margin-right: 4px; }}
  .methodology {{ background: {COLORS["bg_card"]}; border: 1px solid {COLORS["border"]}; border-radius: 8px; padding: 24px; margin: 24px 0; font-size: 0.9em; color: {COLORS["text_muted"]}; }}
  .methodology h3 {{ color: {COLORS["secondary"]}; margin-bottom: 12px; }}
  .methodology p {{ margin: 8px 0; }}
  .disclaimer {{ border-left: 3px solid {COLORS["accent"]}; padding: 12px 16px; margin: 16px 0; font-size: 0.85em; color: {COLORS["text_muted"]}; }}
  .footer {{ text-align: center; padding: 24px; border-top: 1px solid {COLORS["border"]}; color: {COLORS["text_muted"]}; font-size: 0.85em; margin-top: 40px; }}
  @media (max-width: 768px) {{
    .detail-grid {{ grid-template-columns: 1fr; }}
    .container {{ padding: 16px; }}
    .header {{ padding: 24px; }}
  }}
</style>
</head>
<body>

<div class="header">
  <h1><span>Archaeon</span> Extremophile Enzyme Discovery Pipeline</h1>
  <p>Report generated: {timestamp} &nbsp;|&nbsp; Version {version} &nbsp;|&nbsp;
     Sources: {sources} &nbsp;|&nbsp; Total candidates evaluated: {n_total}</p>
</div>

<div class="container">

  <div class="stat-grid">
    <div class="stat-card">
      <div class="stat-value">{len(candidates)}</div>
      <div class="stat-label">Candidates Ranked</div>
    </div>
    <div class="stat-card">
      <div class="stat-value">{max_ias:.1f}</div>
      <div class="stat-label">Top IAS Score</div>
    </div>
    <div class="stat-card">
      <div class="stat-value">{mean_ias:.1f}</div>
      <div class="stat-label">Mean IAS Score</div>
    </div>
    <div class="stat-card">
      <div class="stat-value">{structure_count}</div>
      <div class="stat-label">Structures Predicted</div>
    </div>
    <div class="stat-card">
      <div class="stat-value">{top_family}</div>
      <div class="stat-label">Top Predicted Family</div>
    </div>
    <div class="stat-card">
      <div class="stat-value" style="font-size:1.2em">{top_organism[:20]}</div>
      <div class="stat-label">Top Candidate Organism</div>
    </div>
  </div>

  <div class="disclaimer">
    <strong>Research Tool Disclaimer:</strong> Archaeon is a computational discovery pipeline.
    All scores and predictions are theoretical. Candidates require experimental validation
    before any industrial, clinical, or commercial use. Percent identity to reference enzymes
    does not guarantee functional equivalence.
  </div>

  <div class="section-title">Score Distribution and Analysis</div>
  <div class="charts-grid">
    <div class="chart-card"><canvas id="chart-distribution"  style="height:320px;"></canvas></div>
    <div class="chart-card"><canvas id="chart-top10"  style="height:320px;"></canvas></div>
    <div class="chart-card"><canvas id="chart-scatter"  style="height:320px;"></canvas></div>
  </div>

  <div class="section-title">Ranked Candidates</div>
  <div style="overflow-x:auto;">
    <table class="candidates-table">
      <thead>
        <tr>
          <th>Rank</th><th>Accession</th><th>Description</th><th>Organism</th>
          <th>Biome</th><th>Predicted Family</th><th>IAS</th>
          <th>BLAST Identity</th><th>Length (aa)</th>
        </tr>
      </thead>
      <tbody>
        {table_rows}
      </tbody>
    </table>
  </div>

  {"<div class='section-title'>Predicted Structures (ESMFold)</div>" + structure_viewers if structure_viewers else ""}

  <div class="methodology">
    <h3>Methodology Notes</h3>
    <p><strong>Data Sources:</strong> Protein sequences were retrieved from NCBI Entrez
       (protein database) and MGnify (metagenomic protein database) using targeted biome
       queries for extreme environments including hydrothermal vents, hot springs,
       hypersaline lakes, and acidic environments.</p>
    <p><strong>IAS Formula:</strong> IAS = 0.40 * Thermostability + 0.20 * Quality + 0.40 * BLAST.
       Thermostability component combines aliphatic index (30%), instability index (25%),
       charged residue fraction (25%), proline content (10%), and aromaticity (10%).
       BLAST component uses percent identity to the best match among 8 reference
       thermostable enzyme families.</p>
    <p><strong>Structure Prediction:</strong> ESMFold (Meta AI, Lin et al. 2022) via public
       REST API. pLDDT scores embedded in B-factor column of PDB output. Only sequences
       under 400 aa predicted (truncated otherwise).</p>
    <p><strong>Limitations:</strong> IAS scores are theoretical rankings, not experimental
       thermostability measurements. BLAST identity cutoff for reliable function annotation
       is approximately 40-60% depending on protein family. pLDDT scores are confidence
       estimates for structural coordinates, not thermostability predictions.</p>
  </div>

  <div class="footer">
    Generated by Archaeon v{version} &nbsp;|&nbsp;
    Angie Xiu &nbsp;|&nbsp;
    The Xiu Lab &nbsp;|&nbsp;
    {datetime.now().year} &nbsp;|&nbsp;
    Data from NCBI and MGnify (EBI) &nbsp;|&nbsp;
    Structure prediction by ESMFold (Meta AI)
  </div>

</div>

<script>
function toggleDetail(id) {{
  var el = document.getElementById(id);
  if (el) el.style.display = (el.style.display === 'none' || el.style.display === '') ? 'table-row' : 'none';
}}
{charts_js}

function initViewers() {{
  var pdbScripts = document.querySelectorAll('script[type="text/pdb-data"]');
  pdbScripts.forEach(function(script) {{
    var rank = script.id.replace('pdb-data-', '');
    var element = document.getElementById('viewer-' + rank);
    if (!element) return;
    var pdbData = script.textContent;
    var config = {{ backgroundColor: '#0d1117', antialias: true }};
    var viewer = $3Dmol.createViewer(element, config);
    viewer.addModel(pdbData, 'pdb');
    viewer.setStyle({{}}, {{cartoon: {{colorscheme: {{
      prop: 'b',
      gradient: 'linear',
      min: 0, max: 100,
      colors: ['#FF7D45', '#FFDB13', '#65CBF3', '#0053D6']
    }}}}}});
    viewer.zoomTo();
    viewer.render();
  }});
}}

if (typeof $3Dmol !== 'undefined') {{
  initViewers();
}} else {{
  document.querySelector('script[src*="3Dmol"]').addEventListener('load', initViewers);
}}
</script>

</body>
</html>"""

    output_dir = os.path.dirname(output_path)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html)

    logger.info(f"Report written to {output_path}")
    return output_path
