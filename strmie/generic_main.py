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
    parser.add_argument("--nanopore", action="store_true", help="Fuzzy flank-anchored matching for noisy long reads (Oxford Nanopore), tolerant of errors landing inside the repeat tract itself, not just at read termini. Requires the locus to define flank_upstream/flank_downstream (built-in FMR1/C9orf72 do). Default is exact motif/pattern matching, unaffected by this flag or the --np-* options below.")
    parser.add_argument("--np-max-roi", dest="np_max_roi", type=int, default=3000, help="Max region-of-interest length in bp between flanks, used only with --nanopore (default: 3000; raise for extremely long expansions).")
    parser.add_argument("--np-max-edits", dest="np_max_edits", type=int, default=2, help="Max edit distance allowed for each flank match, used only with --nanopore (default: 2).")
    parser.add_argument("--np-seed-len", dest="np_seed_len", type=int, default=6, help="Seed prefilter length for flank matching, used only with --nanopore (0 disables; default: 6).")
    bestmatch_group = parser.add_mutually_exclusive_group()
    bestmatch_group.add_argument("--np-bestmatch", dest="np_bestmatch", action="store_true", default=True, help="Use regex.BESTMATCH (globally optimal flank alignment) with --nanopore. Default. Can be much slower for large --np-max-roi; see --np-no-bestmatch.")
    bestmatch_group.add_argument("--np-no-bestmatch", dest="np_bestmatch", action="store_false", help="Disable regex.BESTMATCH with --nanopore (first-fit fuzzy match instead of globally optimal). Much faster for large --np-max-roi, at the cost of a possibly slightly less precise flank alignment.")
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
            nanopore=args.nanopore, np_max_roi=args.np_max_roi,
            np_max_edits=args.np_max_edits, np_seed_len=args.np_seed_len,
            np_use_bestmatch=args.np_bestmatch,
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
