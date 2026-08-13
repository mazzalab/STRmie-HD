#!/usr/bin/env python
# coding: utf-8
"""
strmie-repeat: catalog-driven genotyping for tandem repeat loci beyond HTT.

A separate CLI entry point from `strmie` (strmie/main.py), sharing none of
its code paths -- this exists to demonstrate that STRmie-HD's alignment-free
counting approach generalizes to other loci without touching, or risking
any regression in, the validated HTT-specific pipeline.
"""

import argparse
import os

import pandas as pd
from colorama import Fore, Style

from strmie.scripts.generic_repeat import genotype_sample
from strmie.scripts.repeat_catalog import load_catalog, list_builtin_loci
from strmie.scripts.generic_html_generator import create_html


def main():
    parser = argparse.ArgumentParser(
        description=Fore.GREEN + Style.BRIGHT
        + "strmie-repeat: catalog-driven tandem repeat genotyping (non-HTT loci)"
        + Style.RESET_ALL,
    )
    parser.add_argument(
        "-f", "--input", nargs="+", required=True,
        help="One or more .fastq.gz/.fasta.gz files, and/or directories containing them.",
    )
    parser.add_argument(
        "--locus", required=True,
        help=f"Built-in locus name ({', '.join(list_builtin_loci())}), or a path to a custom JSON catalog file.",
    )
    parser.add_argument("-o", "--output", required=True, help="Output directory.")
    parser.add_argument("-i", dest="intorno", type=int, default=5, help="Minimum repeat-unit separation between the two called alleles (default: 5).")
    parser.add_argument("-ti", dest="threshold_instability", type=float, default=False, help="Relative peak height threshold for the Instability Index (default: False, recommended value: 0.2).")
    parser.add_argument("-te", dest="threshold_expansion", type=float, default=False, help="Relative peak height threshold for the Expansion Index (default: False, recommended value: 0.03).")
    args = parser.parse_args()

    catalog = load_catalog(args.locus)
    os.makedirs(args.output, exist_ok=True)

    input_files = []
    for p in args.input:
        if os.path.isdir(p):
            for name in sorted(os.listdir(p)):
                if name.endswith(".fastq.gz") or name.endswith(".fasta.gz"):
                    input_files.append(os.path.join(p, name))
        elif os.path.isfile(p):
            input_files.append(p)
        else:
            raise FileNotFoundError(f"Input path not found: {p}")

    rows = []
    reads_by_sample = {}
    for path in input_files:
        sample = os.path.basename(path)
        print(f"Genotyping {sample} against {catalog['gene']} ({catalog['motif']} repeat)...")
        result = genotype_sample(
            path, catalog, intorno=args.intorno,
            ii_threshold=args.threshold_instability, ei_threshold=args.threshold_expansion,
        )
        reads = result["reads"]
        reads_by_sample[sample] = reads
        rows.append({
            "Sample": sample,
            "Gene": catalog["gene"],
            "Disease": catalog.get("disease", ""),
            "Motif": catalog["motif"],
            "Allele_1": result["allele1"],
            "Allele_2": result["allele2"],
            "N_reads_with_motif": len(reads),
            "Instability_Index": result["ii"],
            "Expansion_Index": result["ei"],
        })
        reads.to_csv(os.path.join(args.output, sample + ".repeat_counts.csv"), index=False)

    df = pd.DataFrame(rows)
    report_path = os.path.join(args.output, "Generic_repeat_report.xlsx")
    df.to_excel(report_path, index=False)
    print(f"Report written to {report_path}")

    create_html(args.output, rows, reads_by_sample, catalog)
    print(f"HTML report written to {os.path.join(args.output, 'report.html')}")


if __name__ == "__main__":
    main()
