#!/usr/bin/env python
# coding: utf-8
"""
Self-contained HTML report for the generic (non-HTT) submodule, in the same
visual/interaction style as strmie/scripts/html_generator.py's canonical
STRmie-HD report (cohort table, per-sample detail panel with a Plotly
histogram, glossary drawer) -- but not a copy of it: this is a separate
template built for a single-locus, single-histogram result, so it omits
whatever doesn't apply to a generic locus rather than showing broken or
misleading sections:
  - No CCG histogram/allele columns (HTT-specific second tract).
  - No LOI-CAA/LOI-CCA/DOI interruption-percentage panel (HTT's specific
    interruption taxonomy; FMR1's own interruption signal is folded into
    the curated repeat count itself, see repeat_catalog.py).
  - No manual-correction/export-to-Index_Calculation panel (that workflow
    is tied to strmie's own raw_counts-folder convention, which this
    submodule does not use).
It does carry over the same somatic Instability/Expansion Index metrics
(computed by generic_indices.py, reusing indices.py unmodified) and a
catalog-driven category badge (Normal/Premutation/Pathogenic/...) derived
from each locus's own literature ranges in repeat_catalog.py, instead of
HTT's fixed CAG clinical cutpoints.
"""

import json
import os

import pandas as pd


def _to_num(v):
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


def _category(value, catalog):
    """Category badge from the locus's own literature ranges in the
    catalog (normal_range / premutation_range / pathogenic_range /
    full_mutation_range), not a hardcoded HTT cutpoint."""
    n = _to_num(value)
    if n is None:
        return ("Needs review", "warn")

    def _in_range(rng):
        if not rng:
            return False
        lo, hi = rng
        if lo is not None and n < lo:
            return False
        if hi is not None and n > hi:
            return False
        return True

    if _in_range(catalog.get("normal_range")):
        return ("Normal", "normal")
    if _in_range(catalog.get("full_mutation_range")) or _in_range(catalog.get("pathogenic_range")):
        return ("Pathogenic", "full")
    if _in_range(catalog.get("premutation_range")):
        return ("Premutation", "reduced")
    return ("Intermediate", "intermediate")


def _hist_from_series(series):
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


def build_report_payload(rows, reads_by_sample, catalog):
    samples = []
    for row in rows:
        sample = row["Sample"]
        allele1_raw, allele2_raw = row.get("Allele_1"), row.get("Allele_2")
        is_warning = isinstance(allele1_raw, str) or isinstance(allele2_raw, str)
        cat_label, cat_class = _category(None if is_warning else allele2_raw, catalog)

        reads = reads_by_sample.get(sample)
        hist = _hist_from_series(reads["Repeat_count"]) if reads is not None else {}
        max_observed = _to_num(reads["Repeat_count"].max()) if reads is not None and len(reads) else None

        samples.append({
            "sample": str(sample),
            "warning": bool(is_warning),
            "allele1": _to_num(allele1_raw),
            "allele2": _to_num(allele2_raw),
            "ii": _to_num(row.get("Instability_Index")),
            "ei": _to_num(row.get("Expansion_Index")),
            "max_observed": max_observed,
            "category": cat_label,
            "category_class": cat_class,
            "hist": hist,
            "n_reads": int(row.get("N_reads_with_motif", 0)),
        })

    return {
        "gene": catalog.get("gene", ""),
        "disease": catalog.get("disease", ""),
        "motif": catalog.get("motif", ""),
        "n_samples": len(samples),
        "samples": samples,
    }


_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>strmie-repeat Report</title>
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
  input[type="text"] {
    border: 1px solid var(--border);
    border-radius: 6px;
    padding: 8px 10px;
    font-size: 14px;
    font-family: inherit;
    width: 260px;
  }
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
  button.icon-btn {
    width: 34px; height: 34px; padding: 0; border-radius: 50%;
    font-weight: 600; display: inline-flex; align-items: center; justify-content: center;
  }
  button:hover { background: var(--accent-light); }
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
  .muted { color: var(--muted); }
  .genotype { font-variant-numeric: tabular-nums; font-weight: 600; }

  #detailContent { display: none; }
  #detailContent.open { display: block; }
  #detailPlaceholder.hidden { display: none; }
  .detail-header { display: flex; justify-content: space-between; align-items: baseline; margin-bottom: 10px; }
  .detail-header h2 { margin: 0; }
  .metric-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(190px, 1fr)); gap: 14px; margin-bottom: 18px; }
  .metric-card { border: 1px solid var(--border); border-radius: 8px; padding: 12px 14px; background: #fafbfc; }
  .metric-card .label { font-size: 12px; color: var(--muted); }
  .metric-card .value { font-size: 22px; font-weight: 700; margin-top: 3px; }
  .metric-card .interp { font-size: 12.5px; color: var(--muted); margin-top: 6px; }
  .chart-box { width: 100%; max-width: 100%; overflow: hidden; border: 1px solid var(--border); border-radius: 8px; padding: 12px; box-sizing: border-box; margin-bottom: 18px; }
  .chart-box h3 { margin: 0 0 8px; font-size: 13.5px; }
  .empty-state { color: var(--muted); font-size: 13.5px; padding: 20px; text-align: center; }
  .banner {
    background: #fff7ed; border: 1px solid #fdba74; color: #9a3412;
    border-radius: 8px; padding: 10px 14px; font-size: 13px; margin-bottom: 14px;
  }

  #helpOverlay { position: fixed; inset: 0; background: rgba(15,23,42,0.35); display: none; z-index: 30; }
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
</style>
</head>
<body>

<header>
  <div>
    <h1>strmie-repeat Report &mdash; <span id="geneName"></span></h1>
    <div class="meta"><span id="diseaseName"></span> &middot; <span id="motifName"></span> repeat &middot; <span id="sampleCount"></span> samples</div>
  </div>
  <div class="header-actions">
    <button class="icon-btn" id="helpBtn" title="Glossary and help">?</button>
  </div>
</header>

<div class="container">

  <div class="banner">
    This report is generated by <code>strmie-repeat</code>, a separate, catalog-driven submodule
    demonstrating that STRmie-HD's core alignment-free counting approach generalizes beyond HTT. It does
    not include interruption-variant tracking (LOI/DOI) or a CCG-style second tract &mdash; those are
    HTT-specific; see the built-in catalog entry for whether this locus has its own curated,
    interruption-aware pattern.
  </div>

  <div class="card">
    <h2>Cohort overview</h2>
    <p class="muted" style="font-size:13px; margin-top:-6px;">Click any row below to see its full detail and repeat-length histogram.</p>
    <div class="toolbar">
      <input type="text" id="searchBox" placeholder="Search sample name...">
      <span class="muted" id="filteredCount" style="font-size:13px;"></span>
    </div>
    <div class="table-wrap">
      <table id="cohortTable">
        <thead>
          <tr>
            <th data-key="sample">Sample</th>
            <th data-key="allele2">Genotype (Allele 1 / 2)</th>
            <th data-key="category_class">Category</th>
            <th data-key="ii">Instability Index (II)</th>
            <th data-key="ei">Expansion Index (EI)</th>
            <th data-key="n_reads">Reads</th>
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
        <button id="closeDetailBtn">Close</button>
      </div>
      <div id="detailWarningBanner"></div>
      <div class="metric-grid" id="metricGrid"></div>
      <div class="chart-box">
        <h3 id="chartTitle"></h3>
        <div id="repeatChart" style="width:100%;height:280px;"></div>
      </div>
    </div>
  </div>

</div>

<div id="helpOverlay"></div>
<div id="helpDrawer">
  <button class="close-drawer" id="closeHelpBtn">Close</button>
  <h2>Glossary</h2>
  <div class="help-item">
    <div class="term">Allele 1 / Allele 2</div>
    <div class="def">The two called repeat lengths (allele sizes), automatically detected as the two most prominent peaks in the per-sample repeat-length histogram (the same peak-calling algorithm used by STRmie-HD's own HTT pipeline).</div>
  </div>
  <div class="help-item">
    <div class="term">Category</div>
    <div class="def">Classification of the longer (Allele 2) call against this locus's own published size ranges (see the catalog entry), not a fixed HTT cutpoint.</div>
  </div>
  <div class="help-item">
    <div class="term">Instability Index (II)</div>
    <div class="def">Measures the asymmetry of the repeat-length distribution around the called Allele 2. Positive values indicate a net bias toward further expansion (somatic mosaicism skewed toward longer repeats); negative values indicate a net bias toward contraction.</div>
  </div>
  <div class="help-item">
    <div class="term">Expansion Index (EI)</div>
    <div class="def">Quantifies the accumulation of signal beyond the Allele 2 peak (reads longer than the called allele). Higher values indicate a greater degree of somatic expansion beyond the primary allele.</div>
  </div>
  <div class="help-item">
    <div class="term">Max repeat observed</div>
    <div class="def">The single longest repeat length seen in any read for this sample, regardless of the called allele peaks.</div>
  </div>
  <div class="help-item">
    <div class="term">Needs review / warning samples</div>
    <div class="def">Automatic peak detection did not find two clear alleles for this sample (e.g. too few reads, or an ambiguous distribution).</div>
  </div>
</div>

<script id="report-data" type="application/json">__REPORT_DATA_JSON__</script>

<script>
(function () {
  "use strict";

  const RAW = JSON.parse(document.getElementById("report-data").textContent);
  const SAMPLES = RAW.samples;

  document.getElementById("geneName").textContent = RAW.gene;
  document.getElementById("diseaseName").textContent = RAW.disease;
  document.getElementById("motifName").textContent = RAW.motif;
  document.getElementById("sampleCount").textContent = RAW.n_samples;

  let sortKey = "sample";
  let sortDir = 1;
  let filterText = "";
  let selectedSample = null;

  function fmt(n, digits) {
    if (n === null || n === undefined) return "warning";
    if (digits === undefined) digits = 2;
    return Number(n).toFixed(digits).replace(/\.00$/, "");
  }

  function genotypeText(s) {
    if (s.warning) return '<span class="muted">warning: peaks not detected</span>';
    return '<span class="genotype">' + fmt(s.allele1, 0) + ' / ' + fmt(s.allele2, 0) + '</span>';
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
        "<td>" + s.n_reads + "</td>";
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

    const trace = { x: xs, y: ys, type: "bar", marker: { color: colors }, hovertemplate: "%{x} repeats: %{y} reads<extra></extra>" };

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
    const config = { responsive: true, displaylogo: false, modeBarButtonsToRemove: ["lasso2d", "select2d", "autoScale2d"] };
    Plotly.react(containerEl, [trace], layout, config);
  }

  function iiInterpretation(ii) {
    if (ii === null) return "Automatic peak detection failed for this sample; index not computed.";
    if (ii > 0.5) return "Net bias toward expansion: somatic mosaicism skewed toward longer repeats.";
    if (ii < -0.5) return "Net bias toward contraction: distribution skewed toward shorter repeats.";
    return "Roughly balanced distribution around the called allele.";
  }

  function eiInterpretation(ei) {
    if (ei === null) return "Automatic peak detection failed for this sample; index not computed.";
    if (ei > 20) return "Substantial accumulation of reads beyond the called allele: high somatic expansion.";
    if (ei > 5) return "Moderate accumulation of reads beyond the called allele.";
    return "Little to no signal beyond the called allele.";
  }

  function metricCard(label, value, interp) {
    return '<div class="metric-card"><div class="label">' + label + '</div>' +
      '<div class="value">' + value + '</div>' +
      '<div class="interp">' + interp + '</div></div>';
  }

  function openDetail(sampleName) {
    const s = SAMPLES.find(x => x.sample === sampleName);
    if (!s) return;
    selectedSample = sampleName;
    renderCohort();

    document.getElementById("detailSampleName").textContent = s.sample;

    const banner = document.getElementById("detailWarningBanner");
    banner.innerHTML = s.warning
      ? '<div class="banner">Automatic peak detection did not find two clear alleles for this sample (' + s.n_reads + ' reads analyzed). Inspect the histogram below.</div>'
      : "";

    const grid = document.getElementById("metricGrid");
    grid.innerHTML =
      metricCard("Genotype (Allele 1 / 2)", s.warning ? "warning" : (fmt(s.allele1, 0) + " / " + fmt(s.allele2, 0)), s.n_reads + " reads analyzed") +
      metricCard("Category", '<span class="badge ' + s.category_class + '">' + s.category + '</span>', "") +
      metricCard("Max repeat observed", s.max_observed === null ? "warning" : fmt(s.max_observed, 0), "Longest single read, any allele") +
      metricCard("Instability Index (II)", s.ii === null ? "warning" : fmt(s.ii), iiInterpretation(s.ii)) +
      metricCard("Expansion Index (EI)", s.ei === null ? "warning" : fmt(s.ei), eiInterpretation(s.ei));

    document.getElementById("chartTitle").textContent = RAW.motif + " repeat distribution";

    document.getElementById("detailPlaceholder").classList.add("hidden");
    document.getElementById("detailContent").classList.add("open");

    renderBarChart(document.getElementById("repeatChart"), s.hist, [s.allele1, s.allele2], RAW.motif + " repeats");

    document.getElementById("detailPanel").scrollIntoView({ behavior: "smooth", block: "nearest" });
  }

  document.getElementById("closeDetailBtn").addEventListener("click", () => {
    selectedSample = null;
    document.getElementById("detailContent").classList.remove("open");
    document.getElementById("detailPlaceholder").classList.remove("hidden");
    renderCohort();
  });

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
})();
</script>

</body>
</html>
"""


def create_html(outdir, rows, reads_by_sample, catalog):
    payload = build_report_payload(rows, reads_by_sample, catalog)
    payload_json = json.dumps(payload).replace("</script>", "<\\/script>")
    html = _TEMPLATE.replace("__REPORT_DATA_JSON__", payload_json)
    with open(os.path.join(outdir, "report.html"), "w") as f:
        f.write(html)
