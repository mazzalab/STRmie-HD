import json
import os
import pandas as pd


# ---------------------------------------------------------------------------
# Data preparation: turn the final report dataframe + per-read dataframe into
# a small, JSON-serializable payload that gets embedded directly in the HTML
# report at generation time. The report is therefore fully self-contained:
# no manual re-uploading of the Excel report or the raw_counts folder is
# needed to view results, browse per-sample histograms, or flag samples for
# manual review.
# ---------------------------------------------------------------------------

def _to_num(v):
    """Best-effort conversion to a plain float, or None for warning/missing values."""
    if v is None:
        return None
    if isinstance(v, str):
        return None
    try:
        if pd.isna(v):
            return None
    except (TypeError, ValueError):
        pass
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _clinical_category(cag2):
    """Clinical CAG categories for the expanded allele, matching the cutpoints
    already used elsewhere in STRmie-HD (<=26 normal; 27-35 intermediate;
    36-39 reduced penetrance; >=40 full penetrance)."""
    n = _to_num(cag2)
    if n is None:
        return ("Needs review", "warn")
    if n <= 26:
        return ("Normal", "normal")
    elif n <= 35:
        return ("Intermediate", "intermediate")
    elif n <= 39:
        return ("Reduced penetrance", "reduced")
    else:
        return ("Full penetrance", "full")


def _hist_from_series(series):
    """Value-count histogram as a {repeat_length: read_count} dict. Only the
    aggregated histogram is embedded (never individual reads), so payload
    size stays tiny regardless of raw sequencing depth."""
    if series is None or len(series) == 0:
        return {}
    counts = series.dropna().value_counts()
    out = {}
    for k, v in counts.items():
        try:
            key = str(int(k))
        except (TypeError, ValueError):
            key = str(k)
        out[key] = int(v)
    return out


def build_report_payload(final_df, raw_data, cutpoint=27):
    samples = []
    for _, row in final_df.iterrows():
        sample = row.get("Sample")
        cag1_raw = row.get("CAG_repeatsPeak_Allele_1")
        cag2_raw = row.get("CAG_repeatsPeak_Allele_2")
        is_warning = isinstance(cag1_raw, str) or isinstance(cag2_raw, str)
        cat_label, cat_class = _clinical_category(None if is_warning else cag2_raw)

        if raw_data is not None and "filename" in raw_data.columns:
            sample_reads = raw_data[raw_data.filename == sample]
        else:
            sample_reads = None

        cag_hist = _hist_from_series(sample_reads["CAG_repeats"]) if sample_reads is not None and "CAG_repeats" in sample_reads.columns else {}
        ccg_hist = _hist_from_series(sample_reads["CCG_repeats"]) if sample_reads is not None and "CCG_repeats" in sample_reads.columns else {}

        samples.append({
            "sample": str(sample),
            "warning": bool(is_warning),
            "cag1": _to_num(cag1_raw),
            "cag2": _to_num(cag2_raw),
            "ccg1": _to_num(row.get("CCG_allele_1")),
            "ccg2": _to_num(row.get("CCG_allele_2")),
            "max_cag": _to_num(row.get("Max_CAG_observed")),
            "ii": _to_num(row.get("Instability_Index")),
            "ei": _to_num(row.get("Expansion_Index")),
            "allele_ratio": _to_num(row.get("Allele_Ratio")),
            "loi_caa": _to_num(row.get("LOI_CAA")),
            "loi_cca": _to_num(row.get("LOI_CCA")),
            "doi": _to_num(row.get("DOI")),
            "category": cat_label,
            "category_class": cat_class,
            "cag_hist": cag_hist,
            "ccg_hist": ccg_hist,
            "n_reads": int(len(sample_reads)) if sample_reads is not None else 0,
        })

    return {"cutpoint": cutpoint, "n_samples": len(samples), "samples": samples}


# ---------------------------------------------------------------------------
# HTML template. Vanilla JS/CSS throughout; the only external dependency is
# SheetJS (xlsx.js), used solely to write a real .xlsx file for the manual
# allele-correction export. Everything else (tables, sorting, search, bar
# charts) is hand-rolled so the report has no other CDN dependency and no
# jQuery/DataTables version-mismatch failure modes.
# ---------------------------------------------------------------------------

_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>STRmie-HD Report</title>
<script src="https://cdn.jsdelivr.net/npm/xlsx@0.18.5/dist/xlsx.full.min.js"></script>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js" charset="utf-8"></script>
<style>
  :root {
    --bg: #f6f7f9;
    --card-bg: #ffffff;
    --border: #e2e5ea;
    --text: #1f2937;
    --muted: #6b7280;
    --accent: #0f6e6e;
    --accent-light: #e6f2f2;
    --normal: #16a34a;
    --intermediate: #ca8a04;
    --reduced: #ea580c;
    --full: #dc2626;
    --warn: #6b7280;
    --flag: #dc2626;
    --radius: 10px;
  }
  * { box-sizing: border-box; }
  body {
    margin: 0;
    font-family: -apple-system, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    background: var(--bg);
    color: var(--text);
    line-height: 1.45;
  }
  header {
    background: var(--card-bg);
    border-bottom: 1px solid var(--border);
    padding: 18px 28px;
    display: flex;
    align-items: center;
    justify-content: space-between;
    position: sticky;
    top: 0;
    z-index: 20;
  }
  header h1 { font-size: 20px; margin: 0; }
  header .meta { color: var(--muted); font-size: 13px; margin-top: 2px; }
  .header-actions { display: flex; gap: 10px; align-items: center; }
  .container { max-width: 1180px; margin: 0 auto; padding: 24px 28px 80px; }
  .card {
    background: var(--card-bg);
    border: 1px solid var(--border);
    border-radius: var(--radius);
    padding: 20px;
    margin-bottom: 20px;
  }
  .card h2 { margin: 0 0 14px; font-size: 16px; }
  .toolbar { display: flex; gap: 12px; align-items: center; margin-bottom: 14px; flex-wrap: wrap; }
  input[type="text"], input[type="number"] {
    border: 1px solid var(--border);
    border-radius: 6px;
    padding: 8px 10px;
    font-size: 14px;
    font-family: inherit;
  }
  #searchBox { width: 260px; }
  button {
    font-family: inherit;
    font-size: 13px;
    border-radius: 6px;
    border: 1px solid var(--border);
    background: var(--card-bg);
    color: var(--text);
    padding: 8px 14px;
    cursor: pointer;
  }
  button.primary { background: var(--accent); border-color: var(--accent); color: white; }
  button.primary:hover { opacity: 0.9; }
  button:hover { background: var(--accent-light); }
  button.icon-btn {
    width: 34px; height: 34px; padding: 0; border-radius: 50%;
    font-weight: 600; display: inline-flex; align-items: center; justify-content: center;
  }
  table { width: 100%; border-collapse: collapse; font-size: 13.5px; }
  thead th {
    text-align: left;
    position: sticky; top: 0;
    background: #fafbfc;
    border-bottom: 1px solid var(--border);
    padding: 10px 12px;
    cursor: pointer;
    user-select: none;
    white-space: nowrap;
  }
  thead th:hover { color: var(--accent); }
  thead th .arrow { color: var(--muted); font-size: 11px; margin-left: 3px; }
  tbody td { padding: 9px 12px; border-bottom: 1px solid var(--border); vertical-align: middle; }
  tbody tr { cursor: pointer; }
  tbody tr:hover { background: var(--accent-light); }
  tbody tr.selected { background: var(--accent-light); }
  .table-wrap { max-height: 520px; overflow: auto; border: 1px solid var(--border); border-radius: 8px; }
  .badge {
    display: inline-block; padding: 3px 10px; border-radius: 999px;
    font-size: 12px; font-weight: 600; color: white; white-space: nowrap;
  }
  .badge.normal { background: var(--normal); }
  .badge.intermediate { background: var(--intermediate); }
  .badge.reduced { background: var(--reduced); }
  .badge.full { background: var(--full); }
  .badge.warn { background: var(--warn); }
  .flag-dot {
    display: inline-block; width: 9px; height: 9px; border-radius: 50%;
    background: var(--flag); margin-right: 5px; vertical-align: middle;
  }
  .muted { color: var(--muted); }
  .genotype { font-variant-numeric: tabular-nums; font-weight: 600; }

  /* Detail panel */
  #detailContent { display: none; }
  #detailContent.open { display: block; }
  #detailPlaceholder.hidden { display: none; }
  .detail-header { display: flex; justify-content: space-between; align-items: baseline; margin-bottom: 10px; }
  .detail-header h2 { margin: 0; }
  .metric-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(190px, 1fr)); gap: 14px; margin-bottom: 18px; }
  .metric-card {
    border: 1px solid var(--border); border-radius: 8px; padding: 12px 14px; background: #fafbfc;
  }
  .metric-card .label { font-size: 12px; color: var(--muted); display: flex; align-items: center; gap: 5px; }
  .metric-card .value { font-size: 22px; font-weight: 700; margin-top: 3px; }
  .metric-card .interp { font-size: 12.5px; color: var(--muted); margin-top: 6px; }
  .info-dot {
    display: inline-flex; align-items: center; justify-content: center;
    width: 15px; height: 15px; border-radius: 50%; background: #d7dde3; color: #45505c;
    font-size: 10px; font-weight: 700; cursor: help;
  }
  .chart-row { display: flex; flex-direction: column; gap: 20px; margin-bottom: 18px; }
  .chart-box { width: 100%; max-width: 100%; overflow: hidden; border: 1px solid var(--border); border-radius: 8px; padding: 12px; box-sizing: border-box; }
  .chart-box h3 { margin: 0 0 8px; font-size: 13.5px; }
  .interruption-row { display: flex; align-items: center; gap: 10px; margin-bottom: 10px; }
  .interruption-label { width: 130px; font-size: 13px; }
  .interruption-track { flex: 1; height: 10px; background: #eceff2; border-radius: 6px; overflow: hidden; }
  .interruption-fill { height: 100%; background: var(--accent); }
  .interruption-fill.flagged { background: var(--flag); }
  .interruption-value { width: 55px; text-align: right; font-size: 13px; font-variant-numeric: tabular-nums; }

  /* Correction list */
  #correctionTable td, #correctionTable th { padding: 8px 10px; }
  #correctionTable input { width: 90px; }
  .empty-state { color: var(--muted); font-size: 13.5px; padding: 20px; text-align: center; }

  /* Help drawer */
  #helpOverlay {
    position: fixed; inset: 0; background: rgba(15,23,42,0.35); display: none; z-index: 30;
  }
  #helpOverlay.open { display: block; }
  #helpDrawer {
    position: fixed; top: 0; right: -420px; width: 400px; max-width: 90vw; height: 100%;
    background: var(--card-bg); box-shadow: -4px 0 24px rgba(0,0,0,0.15);
    transition: right 0.2s ease; z-index: 31; overflow-y: auto; padding: 22px;
  }
  #helpDrawer.open { right: 0; }
  #helpDrawer h2 { margin-top: 0; }
  .help-item { margin-bottom: 18px; }
  .help-item .term { font-weight: 700; font-size: 14px; }
  .help-item .def { font-size: 13.5px; color: var(--muted); margin-top: 4px; }
  .close-drawer { float: right; }
  .banner {
    background: #fff7ed; border: 1px solid #fdba74; color: #9a3412;
    border-radius: 8px; padding: 10px 14px; font-size: 13px; margin-bottom: 14px;
  }
</style>
</head>
<body>

<header>
  <div>
    <h1>STRmie-HD Report</h1>
    <div class="meta"><span id="sampleCount"></span> samples</div>
  </div>
  <div class="header-actions">
    <button class="icon-btn" id="helpBtn" title="Glossary and help">?</button>
  </div>
</header>

<div class="container">

  <div class="card">
    <h2>Cohort overview</h2>
    <p class="muted" style="font-size:13px; margin-top:-6px;">Click any row below to see its full detail, histograms, and to flag it for manual allele correction.</p>
    <div class="toolbar">
      <input type="text" id="searchBox" placeholder="Search sample name...">
      <span class="muted" id="filteredCount" style="font-size:13px;"></span>
    </div>
    <div class="table-wrap">
      <table id="cohortTable">
        <thead>
          <tr>
            <th data-key="sample">Sample</th>
            <th data-key="cag2">Genotype (CAG1 / CAG2)</th>
            <th data-key="category_class">Category</th>
            <th data-key="ii">Instability Index (II)</th>
            <th data-key="ei">Expansion Index (EI)</th>
            <th data-key="allele_ratio">Allele Ratio</th>
            <th>Interruption flags</th>
          </tr>
        </thead>
        <tbody id="cohortBody"></tbody>
      </table>
    </div>
  </div>

  <div class="card" id="detailPanel">
    <div class="empty-state" id="detailPlaceholder">No sample selected. Click a row in the cohort table above to see its detail here.</div>
    <div id="detailContent">
      <div class="detail-header">
        <h2 id="detailSampleName"></h2>
        <div>
          <button id="flagForReviewBtn" class="primary">+ Add to manual review list</button>
          <button id="closeDetailBtn">Close</button>
        </div>
      </div>
      <div id="detailWarningBanner"></div>
      <div class="metric-grid" id="metricGrid"></div>
      <div class="chart-row">
        <div class="chart-box">
          <h3>CAG repeat distribution</h3>
          <div id="cagChart" style="width:100%;height:260px;"></div>
        </div>
        <div class="chart-box">
          <h3>CCG repeat distribution</h3>
          <div id="ccgChart" style="width:100%;height:260px;"></div>
        </div>
      </div>
      <h3 style="font-size:13.5px;">Interruption variants (% of reads)</h3>
      <div id="interruptionBars"></div>
    </div>
  </div>

  <div class="card">
    <h2>Manual allele correction</h2>
    <p class="muted" style="font-size:13px; margin-top:-6px;">
      Add samples here whose automatically detected peaks do not match what you see in the histogram above.
      Fields are pre-filled with the currently detected values so you only need to edit the ones that are wrong.
    </p>
    <div class="toolbar">
      <button id="exportCorrectionsBtn" class="primary">Export corrected alleles (.xlsx)</button>
    </div>
    <div class="banner" style="margin-top:4px;">
      <strong>To apply your corrections:</strong> run STRmie-HD again in Index_Calculation mode, pointing
      <code>-o</code> at <strong>the same output folder</strong> used for this run (it must contain the
      <code>raw_counts</code> subfolder this run created), and <code>-p</code> at the file you just exported:<br>
      <code>strmie --mode Index_Calculation -o &lt;this run's output folder&gt; -p CAG_data_for_recalculating_indices.xlsx</code><br>
      This writes <code>indices_calculation.xlsx</code> with indices recomputed from your corrected allele lengths.
      <code>-f</code>/<code>--input</code> is not needed in this mode.
    </div>
    <table id="correctionTable">
      <thead><tr><th>Sample</th><th>CAG Allele 1</th><th>CAG Allele 2</th><th></th></tr></thead>
      <tbody id="correctionBody"></tbody>
    </table>
    <div class="empty-state" id="correctionEmpty">No samples added yet. Click a sample row in the cohort table above, then click "+ Add to manual review list" in its detail view.</div>
  </div>

</div>

<div id="helpOverlay"></div>
<div id="helpDrawer">
  <button class="close-drawer" id="closeHelpBtn">Close</button>
  <h2>Glossary</h2>
  <div class="help-item">
    <div class="term">CAG_repeatsPeak_Allele_1 / Allele_2</div>
    <div class="def">The two called CAG repeat lengths (allele sizes), automatically detected as the two most prominent peaks in the per-sample CAG-repeat-length histogram.</div>
  </div>
  <div class="help-item">
    <div class="term">Category</div>
    <div class="def">Clinical range for the longer (expanded) allele: Normal &le;26, Intermediate 27&ndash;35, Reduced penetrance 36&ndash;39, Full penetrance &ge;40 CAG repeats.</div>
  </div>
  <div class="help-item">
    <div class="term">Instability Index (II)</div>
    <div class="def">Measures the asymmetry of the CAG distribution around the expanded allele. Positive values indicate a net bias toward further expansion (somatic mosaicism skewed toward longer repeats); negative values indicate a net bias toward contraction. Values near zero indicate a balanced, symmetric distribution.</div>
  </div>
  <div class="help-item">
    <div class="term">Expansion Index (EI)</div>
    <div class="def">Quantifies the accumulation of signal beyond the expanded allele peak (reads longer than the called allele). Higher values indicate a greater degree of somatic expansion beyond the primary allele.</div>
  </div>
  <div class="help-item">
    <div class="term">Allele Ratio</div>
    <div class="def">Ratio of total read signal above the clinical cutpoint (phenotypic zone) to signal at or below the cutpoint (healthy zone). A value of 0 means no reads were observed above the cutpoint; this is expected e.g. for samples with both alleles in the normal range, or can indicate an allele sitting very close to the cutpoint.</div>
  </div>
  <div class="help-item">
    <div class="term">LOI_CAA</div>
    <div class="def">Percentage of analyzed reads showing loss of the canonical CAA interruption in the polyQ tract. LOI-CAA is associated with accelerated disease onset. Flagged (highlighted) above 10%.</div>
  </div>
  <div class="help-item">
    <div class="term">LOI_CCA</div>
    <div class="def">Percentage of analyzed reads showing loss of the canonical CCA motif in the proline-rich (PRD) tract. Flagged above 10%.</div>
  </div>
  <div class="help-item">
    <div class="term">DOI</div>
    <div class="def">Percentage of analyzed reads showing duplication of CAACAG motifs (Duplication of Interruption) in the polyQ tract. Flagged above 10%.</div>
  </div>
  <div class="help-item">
    <div class="term">CCG_allele_1 / Allele_2</div>
    <div class="def">The most common CCG repeat length among reads whose CAG repeat length matches each called CAG allele (within a small tolerance, since the called allele is a histogram peak position and may not exactly equal any single read's CAG count). Shown as "warning" if no reads were found near that allele even with tolerance.</div>
  </div>
  <div class="help-item">
    <div class="term">Max_CAG_observed</div>
    <div class="def">The single longest CAG repeat length seen in any read for this sample, regardless of the called allele peaks. Useful for spotting rare, highly expanded subclones.</div>
  </div>
  <div class="help-item">
    <div class="term">Needs review / warning samples</div>
    <div class="def">Automatic peak detection did not find two clear alleles for this sample (e.g. too few reads, or an ambiguous distribution). Inspect the histogram and use the manual allele correction panel to supply the correct values.</div>
  </div>
  <div class="help-item">
    <div class="term">Index_Calculation mode</div>
    <div class="def">A second STRmie-HD run mode that recomputes II/EI and all other indices from allele lengths you supply yourself, instead of the automatically detected peaks. It reuses the raw per-read counts already saved in this run's <code>raw_counts</code> folder, so it must be pointed (<code>-o</code>) at the same output folder as this run &mdash; it does not need <code>-f</code>/<code>--input</code>. Use the "Manual allele correction" panel below to build and export the required input file.</div>
  </div>
</div>

<script id="report-data" type="application/json">__REPORT_DATA_JSON__</script>

<script>
(function () {
  "use strict";

  const RAW = JSON.parse(document.getElementById("report-data").textContent);
  const SAMPLES = RAW.samples;
  const CUTPOINT = RAW.cutpoint;

  document.getElementById("sampleCount").textContent = RAW.n_samples;

  // ---------------------------------------------------------------------
  // Cohort table: render, sort, search
  // ---------------------------------------------------------------------
  let sortKey = "sample";
  let sortDir = 1;
  let filterText = "";
  let selectedSample = null;
  const corrections = new Map(); // sample -> {cag1, cag2}

  function fmt(n, digits) {
    if (n === null || n === undefined) return "warning";
    if (digits === undefined) digits = 2;
    return Number(n).toFixed(digits).replace(/\.00$/, "");
  }

  function flagDots(s) {
    const flags = [];
    if (s.loi_caa !== null && s.loi_caa > 10) flags.push("LOI-CAA");
    if (s.loi_cca !== null && s.loi_cca > 10) flags.push("LOI-CCA");
    if (s.doi !== null && s.doi > 10) flags.push("DOI");
    if (flags.length === 0) return '<span class="muted">&mdash;</span>';
    return flags.map(f => '<span class="flag-dot"></span>' + f).join("  ");
  }

  function genotypeText(s) {
    if (s.warning) return '<span class="muted">warning: peaks not detected</span>';
    return '<span class="genotype">' + fmt(s.cag1, 0) + ' / ' + fmt(s.cag2, 0) + '</span>';
  }

  function renderCohort() {
    let rows = SAMPLES.filter(s => s.sample.toLowerCase().includes(filterText.toLowerCase()));

    rows.sort((a, b) => {
      let av = a[sortKey], bv = b[sortKey];
      if (av === null || av === undefined) av = -Infinity;
      if (bv === null || bv === undefined) bv = -Infinity;
      if (typeof av === "string") { av = av.toLowerCase(); bv = String(bv).toLowerCase(); }
      if (av < bv) return -1 * sortDir;
      if (av > bv) return 1 * sortDir;
      return 0;
    });

    document.getElementById("filteredCount").textContent = rows.length + " shown";

    const body = document.getElementById("cohortBody");
    body.innerHTML = "";
    rows.forEach(s => {
      const tr = document.createElement("tr");
      if (s.sample === selectedSample) tr.classList.add("selected");
      tr.innerHTML =
        "<td>" + s.sample + "</td>" +
        "<td>" + genotypeText(s) + "</td>" +
        "<td><span class='badge " + s.category_class + "'>" + s.category + "</span></td>" +
        "<td>" + (s.ii === null ? '<span class="muted">warning</span>' : fmt(s.ii)) + "</td>" +
        "<td>" + (s.ei === null ? '<span class="muted">warning</span>' : fmt(s.ei)) + "</td>" +
        "<td>" + (s.allele_ratio === null ? '<span class="muted">warning</span>' : fmt(s.allele_ratio)) + "</td>" +
        "<td>" + flagDots(s) + "</td>";
      tr.addEventListener("click", () => openDetail(s.sample));
      body.appendChild(tr);
    });
  }

  document.querySelectorAll("#cohortTable thead th[data-key]").forEach(th => {
    th.addEventListener("click", () => {
      const key = th.getAttribute("data-key");
      if (sortKey === key) { sortDir *= -1; } else { sortKey = key; sortDir = 1; }
      renderCohort();
    });
  });

  document.getElementById("searchBox").addEventListener("input", (e) => {
    filterText = e.target.value;
    renderCohort();
  });

  // ---------------------------------------------------------------------
  // Interactive Plotly bar chart (hover tooltips, zoom, pan)
  // ---------------------------------------------------------------------
  function renderBarChart(containerEl, histObj, peakValues, axisTitle) {
    const entries = Object.keys(histObj).map(k => [parseFloat(k), histObj[k]]).sort((a, b) => a[0] - b[0]);

    if (entries.length === 0) {
      containerEl.innerHTML = '<div class="muted" style="padding:90px 0; text-align:center; font-size:13px;">No read-level data available</div>';
      return;
    }

    const xs = entries.map(e => e[0]);
    const ys = entries.map(e => e[1]);
    const isPeak = (v) => peakValues && peakValues.some(p => p !== null && Math.round(p) === Math.round(v));
    const colors = xs.map(v => isPeak(v) ? "#dc2626" : "#0f6e6e");

    const trace = {
      x: xs,
      y: ys,
      type: "bar",
      marker: { color: colors },
      hovertemplate: "%{x} repeats: %{y} reads<extra></extra>",
    };

    const annotations = [];
    if (peakValues) {
      peakValues.filter(p => p !== null).forEach((p, i) => {
        annotations.push({
          x: Math.round(p), y: 0, yref: "paper", yshift: -6,
          text: "Allele " + (i + 1), showarrow: false, yanchor: "top",
          font: { size: 10, color: "#dc2626" },
        });
      });
    }

    const layout = {
      margin: { l: 40, r: 10, t: 10, b: annotations.length ? 34 : 26 },
      xaxis: { title: axisTitle || "Repeat length", tickfont: { size: 10 }, fixedrange: false },
      yaxis: { title: "Reads", tickfont: { size: 10 } },
      annotations: annotations,
      font: { family: "-apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif", size: 11 },
      plot_bgcolor: "#fafbfc",
      paper_bgcolor: "#fafbfc",
      bargap: 0.15,
    };

    const config = {
      responsive: true,
      displaylogo: false,
      modeBarButtonsToRemove: ["lasso2d", "select2d", "autoScale2d"],
    };

    Plotly.react(containerEl, [trace], layout, config);
  }

  // ---------------------------------------------------------------------
  // Detail panel
  // ---------------------------------------------------------------------
  function iiInterpretation(ii) {
    if (ii === null) return "Automatic peak detection failed for this sample; index not computed.";
    if (ii > 0.5) return "Net bias toward expansion: somatic mosaicism skewed toward longer repeats.";
    if (ii < -0.5) return "Net bias toward contraction: distribution skewed toward shorter repeats (check for possible PCR under-amplification of the longer allele).";
    return "Roughly balanced distribution around the expanded allele.";
  }

  function eiInterpretation(ei) {
    if (ei === null) return "Automatic peak detection failed for this sample; index not computed.";
    if (ei > 20) return "Substantial accumulation of reads beyond the expanded allele: high somatic expansion.";
    if (ei > 5) return "Moderate accumulation of reads beyond the expanded allele.";
    return "Little to no signal beyond the expanded allele.";
  }

  function ratioInterpretation(r) {
    if (r === null) return "Automatic peak detection failed for this sample; index not computed.";
    if (r === 0) return "No reads observed above the cutpoint (" + CUTPOINT + " CAG), expected for samples with both alleles in the normal range.";
    return "Ratio of phenotypic-zone to healthy-zone signal, split at " + CUTPOINT + " CAG.";
  }

  function metricCard(label, value, interp) {
    return '<div class="metric-card"><div class="label">' + label + '</div>' +
      '<div class="value">' + value + '</div>' +
      '<div class="interp">' + interp + '</div></div>';
  }

  function interruptionRow(label, value) {
    const v = value === null ? 0 : value;
    const flagged = v > 10;
    return '<div class="interruption-row">' +
      '<div class="interruption-label">' + label + '</div>' +
      '<div class="interruption-track"><div class="interruption-fill' + (flagged ? ' flagged' : '') + '" style="width:' + Math.min(v, 100) + '%"></div></div>' +
      '<div class="interruption-value">' + (value === null ? "warning" : v.toFixed(1) + '%') + '</div>' +
      '</div>';
  }

  function openDetail(sampleName) {
    const s = SAMPLES.find(x => x.sample === sampleName);
    if (!s) return;
    selectedSample = sampleName;
    renderCohort();

    document.getElementById("detailSampleName").textContent = s.sample;

    const banner = document.getElementById("detailWarningBanner");
    banner.innerHTML = s.warning
      ? '<div class="banner">Automatic peak detection did not find two clear alleles for this sample (' + s.n_reads + ' reads analyzed). Inspect the histogram below and add this sample to the manual review list.</div>'
      : "";

    const grid = document.getElementById("metricGrid");
    grid.innerHTML =
      metricCard("Genotype (CAG1 / CAG2)", s.warning ? "warning" : (fmt(s.cag1, 0) + " / " + fmt(s.cag2, 0)), s.n_reads + " reads analyzed") +
      metricCard("Category", '<span class="badge ' + s.category_class + '">' + s.category + '</span>', "") +
      metricCard("CCG (Allele 1 / 2)", s.warning ? "warning" : (fmt(s.ccg1, 0) + " / " + fmt(s.ccg2, 0)), "") +
      metricCard("Max CAG observed", s.max_cag === null ? "warning" : fmt(s.max_cag, 0), "Longest single read, any allele") +
      metricCard("Instability Index (II)", s.ii === null ? "warning" : fmt(s.ii), iiInterpretation(s.ii)) +
      metricCard("Expansion Index (EI)", s.ei === null ? "warning" : fmt(s.ei), eiInterpretation(s.ei)) +
      metricCard("Allele Ratio", s.allele_ratio === null ? "warning" : fmt(s.allele_ratio), ratioInterpretation(s.allele_ratio));

    document.getElementById("interruptionBars").innerHTML =
      interruptionRow("LOI-CAA", s.loi_caa) +
      interruptionRow("LOI-CCA", s.loi_cca) +
      interruptionRow("DOI", s.doi);

    // Make the panel visible before rendering the charts: Plotly measures the
    // container's actual pixel width at render time, and a display:none
    // ancestor reports zero width, which produced badly-sized/overflowing charts.
    document.getElementById("detailPlaceholder").classList.add("hidden");
    document.getElementById("detailContent").classList.add("open");

    renderBarChart(document.getElementById("cagChart"), s.cag_hist, [s.cag1, s.cag2], "CAG repeats");
    renderBarChart(document.getElementById("ccgChart"), s.ccg_hist, [s.ccg1, s.ccg2], "CCG repeats");

    document.getElementById("detailPanel").scrollIntoView({ behavior: "smooth", block: "nearest" });
  }

  document.getElementById("closeDetailBtn").addEventListener("click", () => {
    selectedSample = null;
    document.getElementById("detailContent").classList.remove("open");
    document.getElementById("detailPlaceholder").classList.remove("hidden");
    renderCohort();
  });

  // ---------------------------------------------------------------------
  // Manual correction list
  // ---------------------------------------------------------------------
  function renderCorrections() {
    const body = document.getElementById("correctionBody");
    body.innerHTML = "";
    corrections.forEach((vals, sample) => {
      const tr = document.createElement("tr");
      tr.innerHTML =
        "<td>" + sample + "</td>" +
        '<td><input type="number" class="corr-cag1" data-sample="' + sample + '" value="' + (vals.cag1 === null ? "" : vals.cag1) + '"></td>' +
        '<td><input type="number" class="corr-cag2" data-sample="' + sample + '" value="' + (vals.cag2 === null ? "" : vals.cag2) + '"></td>' +
        '<td><button class="corr-remove" data-sample="' + sample + '">Remove</button></td>';
      body.appendChild(tr);
    });
    document.getElementById("correctionEmpty").style.display = corrections.size === 0 ? "block" : "none";

    body.querySelectorAll(".corr-cag1").forEach(inp => inp.addEventListener("input", (e) => {
      const sample = e.target.getAttribute("data-sample");
      corrections.get(sample).cag1 = e.target.value === "" ? null : parseFloat(e.target.value);
    }));
    body.querySelectorAll(".corr-cag2").forEach(inp => inp.addEventListener("input", (e) => {
      const sample = e.target.getAttribute("data-sample");
      corrections.get(sample).cag2 = e.target.value === "" ? null : parseFloat(e.target.value);
    }));
    body.querySelectorAll(".corr-remove").forEach(btn => btn.addEventListener("click", (e) => {
      corrections.delete(e.target.getAttribute("data-sample"));
      renderCorrections();
    }));
  }

  document.getElementById("flagForReviewBtn").addEventListener("click", (e) => {
    if (!selectedSample) return;
    const s = SAMPLES.find(x => x.sample === selectedSample);
    if (!corrections.has(selectedSample)) {
      corrections.set(selectedSample, { cag1: s.cag1, cag2: s.cag2 });
    }
    renderCorrections();
    const btn = e.currentTarget;
    const original = btn.textContent;
    btn.textContent = "✓ Added, see Manual allele correction below";
    btn.disabled = true;
    setTimeout(() => { btn.textContent = original; btn.disabled = false; }, 1800);
  });

  document.getElementById("exportCorrectionsBtn").addEventListener("click", () => {
    if (corrections.size === 0) {
      alert("Add at least one sample to the manual review list first.");
      return;
    }
    const rows = [["Sample", "CAG_Allele_1", "CAG_Allele_2"]];
    corrections.forEach((vals, sample) => rows.push([sample, vals.cag1, vals.cag2]));
    const ws = XLSX.utils.aoa_to_sheet(rows);
    const wb = XLSX.utils.book_new();
    XLSX.utils.book_append_sheet(wb, ws, "Corrections");
    XLSX.writeFile(wb, "CAG_data_for_recalculating_indices.xlsx");
  });

  // ---------------------------------------------------------------------
  // Help drawer
  // ---------------------------------------------------------------------
  document.getElementById("helpBtn").addEventListener("click", () => {
    document.getElementById("helpDrawer").classList.add("open");
    document.getElementById("helpOverlay").classList.add("open");
  });
  function closeHelp() {
    document.getElementById("helpDrawer").classList.remove("open");
    document.getElementById("helpOverlay").classList.remove("open");
  }
  document.getElementById("closeHelpBtn").addEventListener("click", closeHelp);
  document.getElementById("helpOverlay").addEventListener("click", closeHelp);

  renderCohort();
  renderCorrections();
})();
</script>

</body>
</html>
"""


def create_html(outdir, final_df, raw_data, cutpoint=27):
    payload = build_report_payload(final_df, raw_data, cutpoint=cutpoint)
    payload_json = json.dumps(payload).replace("</script>", "<\\/script>")
    html = _TEMPLATE.replace("__REPORT_DATA_JSON__", payload_json)
    with open(outdir + "/report.html", "w") as f:
        f.write(html)


def write_histogram_spreadsheet(outdir, raw_data, out_name="CAG_CCG_histograms.xlsx"):
    """Write the per-sample CAG and CCG repeat-length histograms (the same
    aggregated data embedded in report.html) to a plain spreadsheet, in tidy
    long format: one row per (Sample, Repeat_Type, Repeat_Length) with its
    read count. This makes the distribution data available for downstream
    analysis independent of the interactive report, and independent of the
    much larger per-read raw_counts CSVs."""
    if raw_data is None or "filename" not in raw_data.columns:
        return

    rows = []
    for sample in raw_data["filename"].unique():
        sample_reads = raw_data[raw_data.filename == sample]
        for repeat_type, col in (("CAG", "CAG_repeats"), ("CCG", "CCG_repeats")):
            if col not in sample_reads.columns:
                continue
            hist = _hist_from_series(sample_reads[col])
            for length, count in sorted(hist.items(), key=lambda kv: float(kv[0])):
                rows.append({
                    "Sample": sample,
                    "Repeat_Type": repeat_type,
                    "Repeat_Length": int(float(length)),
                    "Read_Count": count,
                })

    if not rows:
        return

    pd.DataFrame(rows).to_excel(os.path.join(outdir, out_name), index=False)
