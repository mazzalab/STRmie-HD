#!/usr/bin/env python
# coding: utf-8
"""
Per-allele interruption-variant fractions from a STRmie-HD output directory.

STRmie-HD reports interruption-variant fractions (LOI-CAA, LOI-CCA, DOI) per
sample. Its per-read output (`raw_counts/`) records, for every read, both the
CAG repeat length and the interruption status of that same read, so each
interruption call can also be attributed to an allele: every read is assigned
to the nearer of the two called CAG alleles, provided its CAG length lies
within `--window` repeat units of that allele (reads further from both
alleles, e.g. rare somatic expansions or PCR artefacts, are left unassigned),
and the fraction of reads carrying each variant is computed among the reads
of each allele. A variant is attributed to an allele when it is present in at
least `--min-pct` percent of that allele's reads.

Alleles are taken from `Final_report.xlsx` in the output directory, or from a
curated table (`-p`, same format as for --mode Index_Calculation) when the
automatic calls were manually adjusted.

The same computation runs automatically when strmie is launched with
--allele_interruptions (see strmie/main.py).

Usage:
    python -m strmie.scripts.allele_interruptions -o /path/to/output_dir
    python -m strmie.scripts.allele_interruptions -o /path/to/output_dir -p curated.xlsx
"""
import argparse
import glob
import os

import pandas as pd

VARIANTS = (("LOI-CAA", "LOI_CAA"), ("LOI-CCA", "LOI_CCA"), ("DOI", "DOI"))

# Allele column names: Final_report.xlsx uses the first pair, the curated table
# used by --mode Index_Calculation (exported from the HTML report) the second.
ALLELE_COLUMNS = (("CAG_repeatsPeak_Allele_1", "CAG_repeatsPeak_Allele_2"),
                  ("CAG_Allele_1", "CAG_Allele_2"))


def _pct(reads, col):
    if len(reads) == 0:
        return None
    return round(100 * reads[col].astype(str).str.lower().eq("true").mean(), 2)


def allele_fractions(reads, cag1, cag2, window=3, min_pct=50.0):
    """Per-allele interruption fractions for one sample's per-read table."""
    cag1, cag2 = sorted((int(cag1), int(cag2)))
    single = cag1 == cag2
    d1 = (reads.CAG_repeats - cag1).abs()
    d2 = (reads.CAG_repeats - cag2).abs()
    a1 = reads[(d1 <= window) & ((d1 <= d2) | single)]
    a2 = reads.iloc[0:0] if single else reads[(d2 < d1) & (d2 <= window)]

    row = {"CAG_allele_1": cag1, "CAG_allele_2": None if single else cag2,
           "Reads_total": len(reads), "Reads_allele_1": len(a1), "Reads_allele_2": None if single else len(a2)}
    carriers = []
    for label, col in VARIANTS:
        row[f"{label}_%_sample"] = _pct(reads, col)
        p1, p2 = _pct(a1, col), (None if single else _pct(a2, col))
        row[f"{label}_%_allele_1"], row[f"{label}_%_allele_2"] = p1, p2
        on = [name for name, p in (("allele 1", p1), ("allele 2", p2)) if p is not None and p >= min_pct]
        if on:
            carriers.append(f"{label}: {' + '.join(on)}")
    row["Carrying_allele"] = "; ".join(carriers) if carriers else "none"
    return row


def _allele_columns(report):
    for c1, c2 in ALLELE_COLUMNS:
        if c1 in report.columns and c2 in report.columns:
            return c1, c2
    expected = " or ".join(f"{c1}/{c2}" for c1, c2 in ALLELE_COLUMNS)
    raise ValueError(f"Allele table must contain a Sample column and {expected} columns")


def compute_allele_interruptions(report, counts, window=3, min_pct=50.0):
    """Per-allele interruption fractions for every sample of `report`.

    report: one row per sample with a Sample column and the two CAG alleles
            (Final_report.xlsx or an Index_Calculation curated table).
    counts: per-read table with CAG_repeats, LOI_CAA, LOI_CCA, DOI and
            filename columns (the content of raw_counts/).
    Returns a DataFrame with one row per sample."""
    col1, col2 = _allele_columns(report)
    rows = []
    for _, r in report.iterrows():
        sample = r["Sample"]
        cag1, cag2 = r[col1], r[col2]
        reads = counts[counts.filename == sample]
        if isinstance(cag1, str) or isinstance(cag2, str) or pd.isna(cag1) or pd.isna(cag2):
            rows.append({"Sample": sample, "Carrying_allele": "alleles not called"})
            continue
        rows.append({"Sample": sample, **allele_fractions(reads, cag1, cag2, window, min_pct)})
    return pd.DataFrame(rows)


def _read_counts(raw_dir):
    files = sorted(glob.glob(os.path.join(raw_dir, "*.csv")) + glob.glob(os.path.join(raw_dir, "*.csv.gz")))
    if not files:
        raise SystemExit(f"No per-read files found in {raw_dir}")
    cols = ["CAG_repeats", "LOI_CAA", "LOI_CCA", "DOI", "filename"]
    return pd.concat([pd.read_csv(f, usecols=cols) for f in files], ignore_index=True)


def main():
    ap = argparse.ArgumentParser(description="Per-allele interruption-variant fractions from STRmie-HD per-read output.")
    ap.add_argument("-o", "--output", required=True, help="STRmie-HD output directory (contains raw_counts/ and Final_report.xlsx)")
    ap.add_argument("-p", "--peaks", help="Curated allele table (default: <output>/Final_report.xlsx)")
    ap.add_argument("-w", "--window", type=int, default=3, help="Max distance (CAG units) from a called allele for a read to be assigned to it [3]")
    ap.add_argument("--min-pct", type=float, default=50.0, help="Min %% of an allele's reads carrying a variant to attribute it to that allele [50]")
    ap.add_argument("--out", help="Output Excel file (default: <output>/allele_interruptions.xlsx)")
    args = ap.parse_args()

    report = pd.read_excel(args.peaks or os.path.join(args.output, "Final_report.xlsx"))
    counts = _read_counts(os.path.join(args.output, "raw_counts"))
    try:
        result = compute_allele_interruptions(report, counts, args.window, args.min_pct)
    except ValueError as e:
        raise SystemExit(str(e))

    out = args.out or os.path.join(args.output, "allele_interruptions.xlsx")
    result.to_excel(out, index=False)
    print(f"Per-allele interruption fractions for {len(result)} samples written to {out}")


if __name__ == "__main__":
    main()
