#!/usr/bin/env python
# coding: utf-8

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.signal as signal
import math
import os
import shutil
import argparse
from concurrent.futures import ProcessPoolExecutor
from scipy import signal
from scipy.signal import find_peaks


from scipy import stats

from strmie.scripts.utility import *
from strmie.scripts.html_generator import *
from strmie.scripts.indices import *
from strmie.scripts.pattern import *
from strmie.scripts.peaks import *
from strmie.scripts.bam_extract import (
    extract_sample_to_fastq,
    resolve_bam_cram_paths,
)


from colorama import Fore, Style


def _process_one_input_file(task):
    """Worker for parallel per-sample processing in the Complete_Pipeline
    loop below. Kept at module scope so it can be pickled to worker
    processes. Runs exactly the same calcola_counts_and_loi/
    calcola_counts_and_loi_nanopore call the sequential path always used;
    the only difference is SystemExit (raised by those functions on a
    sample with zero CAG repeats found) is converted to a regular
    exception, since SystemExit does not propagate cleanly out of a
    ProcessPoolExecutor worker -- the net effect (the whole run aborts
    with the same message) is unchanged, whether run in parallel or not."""
    path_file, name, nanopore_mode, nanopore_kwargs = task
    try:
        if nanopore_mode:
            df = calcola_counts_and_loi_nanopore(path_file, name=name, **nanopore_kwargs)
        else:
            df = calcola_counts_and_loi(path_file, name)
    except SystemExit as e:
        raise RuntimeError(str(e)) from None
    return name, df


def main():

    parser = argparse.ArgumentParser(description=Fore.GREEN + Style.BRIGHT + '''
    ################################################################################################################
    ##                                                                          .::::-. .=-: ..                   ##
    ##                                                                       .:--:.  .. -=-:.                     ##
    ##                                                           ..:::::::----:::--::::-=+-                       ##
    ##                                                      .::-::::.    .-=-::::-:::---                          ##
    ##                                                   .:--:.    .:    .-=:::                                   ##
    ##                                                 .--:..::.     ::  ==-                                      ##
    ##                                 .....:::.:::::.:=-:.    ....    :.++.                                      ##
    ##                          ..:---------:::--::::---==---------==-:.:*%=                                      ##
    ##                       .-==---::-:       :-.  -===:.     ...::-==++**-                                      ##
    ##                    .-=--:..    .::       := :=+:                                                           ##
    ##                    =--:--.       .::.     -===:                                                            ##
    ##                    -:   .--:       .::     ===                                                             ##
    ##                            :--.      .-:  =++.                                                             ##
    ##                    +:        .--.      :=+++.                                                              ##
    ##                    =#+:        .--:.:-=+++-.                                                               ##
    ##                    C - A - G - C - A - A - C - A - G - C - C - G - C - C - A - C - C - G                   ##
    ##                    |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |                   ##
    ##                    G - T - C - G - T - T - G - T - C - G - G - C - G - G - T - G - G - C                   ##
    ##                                                                                                            ##
    ##               "STRmie-HD - Analyze Huntington's Disease genetic markers from sequencing data"              ##
    ##                                                                                                            ##
    ################################################################################################################
    ''' + Style.RESET_ALL + Fore.MAGENTA + "usage: strmie.py --mode {Complete_Pipeline,Index_Calculation} -f /path/to/raw_reads_file.fastq.gz -o /path/directory_output/ [arguments]" + Style.RESET_ALL , add_help=True, formatter_class=argparse.RawDescriptionHelpFormatter)

    # Definizione della modalità
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--mode", choices=["Complete_Pipeline", "Index_Calculation"], help=Fore.CYAN + "Choose the mode: 'Complete_Pipeline' for full analysis, 'Index_Calculation' for index-based analysis" + Style.RESET_ALL)

    # Parametri per la modalità 'Complete_Pipeline'
    default_group = parser.add_argument_group("Complete Pipeline Arguments")
    default_group.add_argument('-f', '--input', action="store", type=str, nargs='+', required=False, help=Fore.CYAN + "Specify one or more input paths: a directory containing raw reads in .fastq.gz/.fasta.gz format, and/or individual file paths (space-separated). Use this to process a whole folder, a single file, or an explicit list of files. Exactly one of -f/--bam/--cram is required for Complete_Pipeline mode; not used in Index_Calculation mode" + Style.RESET_ALL)
    default_group.add_argument('--bam', action="store", type=str, nargs='+', required=False, help=Fore.CYAN + "Specify one or more indexed BAM files (or directories containing .bam files) instead of FASTQ input. Reads near the HTT CAG/CCG locus are extracted internally (including mates that map away from the locus) and converted to FASTQ before running the normal pipeline unchanged. An exact alternative to -f/--input, not combinable with it." + Style.RESET_ALL)
    default_group.add_argument('--cram', action="store", type=str, nargs='+', required=False, help=Fore.CYAN + "Specify one or more indexed CRAM files (or directories containing .cram files) instead of FASTQ input, same extraction behavior as --bam. Requires --reference." + Style.RESET_ALL)
    default_group.add_argument('--reference', action="store", type=str, required=False, help=Fore.CYAN + "Reference FASTA used to decode --cram input (required with --cram; ignored otherwise)." + Style.RESET_ALL)
    default_group.add_argument('--locus-build', dest='locus_build', choices=["auto", "grch38", "grch37"], default="auto", help=Fore.CYAN + "Reference build for --bam/--cram locus extraction (default: auto-detect from the chr4/4 contig length in the file header)." + Style.RESET_ALL)
    default_group.add_argument('--bam-region', dest='bam_region', action="store", type=str, default=None, help=Fore.CYAN + "Override the --bam/--cram extraction region manually as chrom:start-end (1-based, inclusive), bypassing build auto-detection. Useful for non-HTT loci or non-standard references." + Style.RESET_ALL)
    default_group.add_argument('-o', '--output', action="store", type=str, required=True, help=Fore.CYAN + "Specify the directory path where the output files will be saved (Required)" + Style.RESET_ALL)
    default_group.add_argument('-bc', '--cutpoint_based', action="store_true", default=False, help=Fore.CYAN + "Enable detection of CAG repeat numbers based on identifying the highest peaks in each of the two histograms formed by the cutpoint parameter" + Style.RESET_ALL)
    default_group.add_argument('-a', dest='amp', action='store', type=int, nargs='+', default=[5,6,7,8,9,10], help=Fore.CYAN + 'Specify amplitude values for the expected peak widths in the data (default: [5,6,7,8,9,10])' + Style.RESET_ALL)
    default_group.add_argument('-i', dest='interv', action='store', type=int, default=5, help=Fore.CYAN + 'Minimum CAG-unit separation between the two allele peaks; closer local maxima are treated as noise/stutter of the same allele (default: 5)' + Style.RESET_ALL)
    default_group.add_argument('-ti', dest='threshold_instability', action='store', type=float, default=False, help=Fore.CYAN + 'Set the relative peak height threshold for the instability index (Default False, recommended value: 0.2)' + Style.RESET_ALL)
    default_group.add_argument('-te', dest='threshold_expansion', action='store', type=float, default=False, help=Fore.CYAN + 'Set the relative peak height threshold for the Expansion index (Default False, recommended value: 0.03)' + Style.RESET_ALL)
    default_group.add_argument('-c', dest='cutpoint', action='store', type=int, default=27, help=Fore.CYAN + 'Set the cutpoint to divide the histogram into two sections (default: 27), also used to calculate AlleleRatio ' + Style.RESET_ALL)
    default_group.add_argument('-m', dest='min', action='store', type=int, default=7, help=Fore.CYAN + 'min CAG repeats (default: 7)' + Style.RESET_ALL)
    default_group.add_argument('-j', '--jobs', dest='jobs', action='store', type=int, default=None, help=Fore.CYAN + 'Number of samples to process in parallel during raw-read parsing (default: automatic, up to 4 or the number of input files, whichever is smaller). Set to 1 to process samples sequentially, one at a time. Output is identical regardless of this setting; it only affects wall-clock runtime.' + Style.RESET_ALL)
    default_group.add_argument('--cag_graph', dest='cag', action='store_true', default=False, help=Fore.CYAN + 'Enable to save graphs of CAG trinucleotide repeat distributions' + Style.RESET_ALL)
    default_group.add_argument('--ccg_graph', dest='ccg', action='store_true', default=False, help=Fore.CYAN + 'Enable to save graphs of CCG trinucleotide repeat distributions' + Style.RESET_ALL)
        ##adding nanopore
    # --- Nanopore mode (switch) ---
    default_group.add_argument(
        "--nanopore",
        action="store_true",
        default=False,
        help=Fore.CYAN + "Use Nanopore ROI-based parsing (fuzzy flanks) instead of the default exact-regex method" + Style.RESET_ALL
    )
    ## end adding nanopore
    default_group.add_argument('--cwt', dest='cwt_finder', action='store_true', default=False, help=Fore.CYAN + 'Enable wavelet-based peak detection as an alternative to histogram-based detection' + Style.RESET_ALL)
    default_group.add_argument('--merge_paired_end', action='store_true', default=False, help=Fore.CYAN + "Detect R1/R2 paired-end file pairs in the input (by filename, e.g. *_R1/*_R2 or *_1/*_2) and merge each pair into a single sample with PEAR before processing. Files not part of a detected pair are processed individually as usual. Off by default (each input file is its own sample). Requires the 'pear' executable on PATH" + Style.RESET_ALL)
    default_group.add_argument('--pe_min_overlap', action='store', type=int, default=10, help=Fore.CYAN + 'Minimum overlap (bp) required by PEAR to merge a read pair, used only with --merge_paired_end (default: 10)' + Style.RESET_ALL)

    # Parametri per la modalità 'Index_Calculation'
    indices_group = parser.add_argument_group("Index Calculation Arguments")
    indices_group.add_argument('-p', '--path', action="store", type=str, required=False, help=Fore.CYAN + "Specify the file (in excel format) with four columns detail: Sample, CAG_Allele_1, CAG_Allele_2" + Style.RESET_ALL)

    ##adding nanopore
    # --- Nanopore ROI parsing arguments (only used with --nanopore) ---
    np_group = parser.add_argument_group("Nanopore ROI parsing arguments (used only with --nanopore)")
    np_group.add_argument("--np-max-roi", type=int, default=300, help="Max ROI length (default: 300)")
    np_group.add_argument("--np-max-edits", type=int, default=2, help="Max edits for both flanks")
    np_group.add_argument("--np-max-edits-left", type=int, default=2, help="Override max edits for upstream flank (default: 2)")
    np_group.add_argument("--np-max-edits-right", type=int, default=3, help="Override max edits for downstream flank (default: 3)")

    np_group.add_argument("--np-seed-len", type=int, default=0, help="Seed prefilter length (0 disables). Suggested 5-7.")
    np_group.add_argument("--np-bestmatch", dest="np_bestmatch", action="store_true", default=True,
                          help="Use regex.BESTMATCH (default: True)")
    np_group.add_argument("--np-no-bestmatch", dest="np_bestmatch", action="store_false",
                          help="Disable regex.BESTMATCH (often faster)")


    np_group.add_argument("--np-min-read-len", type=int, default=50, help="Min read length (default: 50)")

    np_group.add_argument("--np-min-cag-pct", type=float, default=0.7,
                          help="Discard reads if the fraction of in-frame CAG triplets in the selected region is below this threshold (default: 0.70). Set 0 to disable.")
    np_group.add_argument("--np-cag-pct-scope", choices=["roi", "cag_region"], default="cag_region",
                          help="Region used for the CAG-fraction filter: 'roi' = entire ROI; 'cag_region' = ROI prefix before the LOI/DOI motif block.")
    np_group.add_argument("--np-allow-caa", action="store_true",
                          help="When computing the fraction, count CAA as acceptable alongside CAG (i.e., CAG or CAA are considered 'good' triplets).")

    ## end adding nanopore
    
    args = parser.parse_args()

    # Logica per verificare la compatibilità dei parametri in base alla modalità
    n_input_modes = sum(bool(x) for x in (args.input, args.bam, args.cram))
    if args.mode == "Complete_Pipeline" and n_input_modes == 0:
        parser.error("In 'Complete_Pipeline' mode, you must provide exactly one of '-f/--input', '--bam', or '--cram'.")
    if args.mode == "Complete_Pipeline" and n_input_modes > 1:
        parser.error("'-f/--input', '--bam', and '--cram' are mutually exclusive; provide exactly one.")
    if args.mode == "Complete_Pipeline" and args.cram and not args.reference:
        parser.error("--cram requires --reference (a FASTA matching the CRAM's alignment reference).")
    if args.mode == "Index_Calculation" and not args.path:
        parser.error("In 'Index_Calculation' mode, you must provide the '-p' or '--path' parameter with a excel file (.xlsx) containing four columns: Sample, CAG_Allele_1, CAG_Allele_2")

    # Assegnazione delle variabili
    if args.mode == "Complete_Pipeline":
        path = args.output + "/"
        os.makedirs(path, exist_ok=True)

        if args.input:
            input_raw_reads = ", ".join(args.input)
            input_files = resolve_input_paths(args.input)

            if args.merge_paired_end:
                pairs, singles = detect_paired_end_pairs(input_files)
                if pairs:
                    merged_paths, pear_stats = merge_paired_end_with_pear(pairs, path, min_overlap=args.pe_min_overlap)
                    input_files = singles + list(merged_paths.values())
                    pd.DataFrame(pear_stats).to_excel(path + "pear_merge_stats.xlsx", index=False)
                    print("Merged " + str(len(pairs)) + " paired-end sample(s) with PEAR: " + ", ".join(p for p, _, _ in pairs))

        else:
            # --bam or --cram: extract HTT-locus reads to FASTQ, then fall through
            # to the exact same downstream pipeline used for -f/--input.
            is_cram = bool(args.cram)
            ext = ".cram" if is_cram else ".bam"
            raw_paths = args.cram if is_cram else args.bam
            input_raw_reads = ", ".join(raw_paths)
            bam_cram_files = resolve_bam_cram_paths(raw_paths, ext)

            region_override = None
            if args.bam_region:
                chrom_part, coords = args.bam_region.split(":")
                start_1based, end_1based = coords.split("-")
                region_override = (chrom_part, int(start_1based) - 1, int(end_1based))

            extract_dir = os.path.join(path, "extracted_fastq")
            os.makedirs(extract_dir, exist_ok=True)

            paired_files = []
            single_by_sample = {}
            for src in bam_cram_files:
                sample_name = os.path.basename(src)[: -len(ext)]
                result = extract_sample_to_fastq(
                    src, sample_name, extract_dir,
                    reference=args.reference if is_cram else None,
                    locus_build=args.locus_build,
                    region_override=region_override,
                )
                if result["paired"]:
                    paired_files.extend(result["paired"])
                if result["single"]:
                    single_by_sample[sample_name] = result["single"]
                parts = []
                if result["paired"]:
                    parts.append("paired")
                if result["single"]:
                    parts.append("single-end/rescued")
                print(f"Extracted {sample_name}: " + " + ".join(parts))

            # BAM/CRAM-derived paired reads have no other route to get merged
            # (the user has no pre-alignment fastqs to run PEAR on themselves),
            # so always merge R1/R2 pairs produced by extraction, independent
            # of --merge_paired_end (which only governs the -f/--input path).
            pairs, _ = detect_paired_end_pairs(paired_files)
            merged_paths = {}
            if pairs:
                merged_paths, pear_stats = merge_paired_end_with_pear(pairs, path, min_overlap=args.pe_min_overlap)
                pd.DataFrame(pear_stats).to_excel(path + "pear_merge_stats.xlsx", index=False)
                print("Merged " + str(len(pairs)) + " paired-end sample(s) extracted from " + ext + ": " + ", ".join(p for p, _, _ in pairs))

            # Fold single-end/rescued content into the same sample's final
            # fastq rather than treating it as a separate sample: a read
            # rescued purely by sequence content (no position-based anchor)
            # commonly has no mate that was independently rescued too, so it
            # never enters the paired/PEAR path at all -- dropping it here
            # would throw away most of what that rescue mechanism finds.
            input_files = []
            handled_samples = set()
            for sample_name, merged_path in merged_paths.items():
                handled_samples.add(sample_name)
                single_path = single_by_sample.get(sample_name)
                if single_path:
                    with open(merged_path, "ab") as out_f, open(single_path, "rb") as in_f:
                        shutil.copyfileobj(in_f, out_f)
                    print(f"Folded rescued single-end reads into {sample_name}'s merged fastq")
                input_files.append(merged_path)
            for sample_name, single_path in single_by_sample.items():
                if sample_name not in handled_samples:
                    input_files.append(single_path)

        ampiezza = args.amp
        intorno = args.interv
        cag_graph = args.cag
        ccg_andWarning_graph = args.ccg

        ## adding nanopore
        nanopore_mode = args.nanopore
        ## end adding nanopore
        
        cwt =args.cwt_finder
        cutpoint = args.cutpoint
        infMin = args.min
        biological_cutpoint = args.cutpoint_based
        ii_threshold = args.threshold_instability
        ei_threshold = args.threshold_expansion
    elif args.mode == "Index_Calculation":
        path = args.output + "/"
        os.makedirs(path, exist_ok=True)
        index_path = args.path
        cutpoint = args.cutpoint
        ii_threshold = args.threshold_instability
        ei_threshold = args.threshold_expansion
    ####

    def report_to_excel(data,campioni,output,cwt):
        instInd=[]
        expInd=[]
        cag=[]
        iiHeight=[]
        cag_max_alleles_1=[]
        cag_max_alleles_2=[]
        percentage_caa=[]
        percentage_cca=[]
        percentage_doi=[]
        observed_maxCAG=[]
        histogramRatio=[]

        for c in campioni:
            data_campione=data[data.filename==c]
        
            ## modification new part
            if biological_cutpoint:
                cag_max_1, cag_max_2 = fine_maxPeak_hist_generated_bycutPoint(data_campione,cutpoint)
            elif cwt: 
                cag_max_1, cag_max_2 = find_peaks_two_alleles(data_campione,ampiezza=ampiezza,intorno=intorno)  # picchi massimi nei due alleli binomial distribution
            else:
                cag_max_1, cag_max_2 = cag_peaks(data_campione, colonna="CAG_repeats",intorno=5)
            ### fine modification new part

            cag_max_alleles_1.append(cag_max_1)
            cag_max_alleles_2.append(cag_max_2)
        

            if type(cag_max_1)==str: ### la stringa warning, altrimenti non è una stringa
                instInd.append("warning")
                expInd.append("warning")
                cag.append("warning")
                iiHeight.append("warning")
                percentage_cca.append("warning")
                percentage_caa.append("warning")
                percentage_doi.append("warning")
                observed_maxCAG.append("warning")
                histogramRatio.append("warning")

            else:
                df_distrib=create_df_distribution(data_campione)
                observed_maxCAG.append(df_distrib["CAG_repeat"].max())
                ii=instabilityIndex(df_distrib,cag_max_1,cag_max_2,pcrFiltering=ii_threshold)  # instability Index , ti
                ei=expansionIndex(df_distrib,cag_max_1,cag_max_2,pcrFiltering=ei_threshold)    # expansion index ,te
                instInd.append(ii)
                expInd.append(ei)
                histogramRatio.append(histogramRatioIndex(df_distrib,cutpoint))
        
                # Data filtering SQUITIERI
                #df_loi = data_campione[data_campione.CAG_repeats >= cag_max_2]
                #df_loi = data_campione.copy()
            
                # Count LOI values
                loi_counts_caa = data_campione.LOI_CAA.value_counts()
                loi_counts_cca = data_campione.LOI_CCA.value_counts()
                doi_counts = data_campione.DOI.value_counts() # new


                # Check for True and False in the index
                loi_caa = loi_counts_caa.get(True, 0)  # Get count of True, default to 0 if not present
                nonLoi_caa = loi_counts_caa.get(False, 0)  # Get count of False, default to 0 if not present

                loi_cca = loi_counts_cca.get(True, 0)  # Get count of True, default to 0 if not present
                nonLoi_cca = loi_counts_cca.get(False, 0)  # Get count of False, default to 0 if not present

                ## DOI
                doi_get = doi_counts.get(True, 0)  # Get count of True, default to 0 if not present
                non_doi = doi_counts.get(False, 0)  # Get count of False, default to 0 if not present

                ##### PRINTING same stats of LOI on screen
                #reads_in_peak = data_campione[data_campione["CAG_repeats"] == cag_max_2]
                #reads_in_peak_caa = reads_in_peak[reads_in_peak["LOI_CAA"] == True]
                #perc_reads_in_peak_caa = (len(reads_in_peak_caa) / len(reads_in_peak) * 100) if len(reads_in_peak) > 0 else 0
                #print(f"{c} - Percentage of LOI_CAA in the peak ({cag_max_2}): {perc_reads_in_peak_caa:.2f}%")

                #reads_in_peak_allele1 = data_campione[data_campione["CAG_repeats"] == cag_max_1]
                #reads_in_peak_doi = reads_in_peak_allele1[reads_in_peak_allele1["DOI"] == True]
                #perc_reads_in_peak_doi = (len(reads_in_peak_doi) / len(reads_in_peak_allele1) * 100) if len(reads_in_peak_allele1) > 0 else 0
                #print(f"{c} - Percentage of DOI in the peak ({cag_max_1}): {perc_reads_in_peak_doi:.2f}%")
                ##########################################


                # Calculate total counts
                total_caa = loi_caa + nonLoi_caa
                total_cca = loi_cca + nonLoi_cca
                total_doi = doi_get + non_doi

                # Calculate percentages
                percentage_caa.append( (loi_caa / total_caa) * 100 if loi_caa > 0 else 0)
                percentage_cca.append((loi_cca / total_cca) * 100 if loi_cca > 0 else 0)
                percentage_doi.append((doi_get / total_doi) * 100 if doi_get > 0 else 0)


        
        df=pd.DataFrame()
        df["Sample"]=campioni
        df["Instability_Index"]=instInd
        df["Expansion_Index"]=expInd
        df["Allele_Ratio"]=histogramRatio
        df["LOI_CAA"]=percentage_caa   
        df["LOI_CCA"]=percentage_cca
        df["DOI"]=percentage_doi
        df["CAG_repeatsPeak_Allele_1"]=cag_max_alleles_1
        df["CAG_repeatsPeak_Allele_2"]=cag_max_alleles_2
        df["Max_CAG_observed"]=observed_maxCAG

        df.to_excel(str(output),index=False)
    
        return df





    def force_findingPeaks(data,campioni,output):
        instInd=[]
        expInd=[]
        cag=[]
        iiHeight=[]
        cag_max_alleles_1=[]
        cag_max_alleles_2=[]
        percentage_caa=[]
        percentage_cca=[]
        percentage_doi=[]
        observed_maxCAG=[]
        histogramRatio=[]

        for c in campioni:
            data_campione=data[data.filename==c]
        
            cag_max_1, cag_max_2 = force_search(data_campione,intorno=intorno)
            cag_max_alleles_1.append(cag_max_1)
            cag_max_alleles_2.append(cag_max_2)
        
            df_distrib=create_df_distribution(data_campione)
            observed_maxCAG.append(df_distrib["CAG_repeat"].max())

            ii=instabilityIndex(df_distrib,cag_max_1,cag_max_2,pcrFiltering=ii_threshold)  # instability Index
            ei=expansionIndex(df_distrib,cag_max_1,cag_max_2,pcrFiltering=ei_threshold)    # expansion index

            histogramRatio.append(histogramRatioIndex(df_distrib,cutpoint))

            instInd.append(ii)
            expInd.append(ei)

                    # Count LOI values
            loi_counts_caa = data_campione.LOI_CAA.value_counts()
            loi_counts_cca = data_campione.LOI_CCA.value_counts()
            doi_counts = data_campione.DOI.value_counts() # new


            # Check for True and False in the index
            loi_caa = loi_counts_caa.get(True, 0)  # Get count of True, default to 0 if not present
            nonLoi_caa = loi_counts_caa.get(False, 0)  # Get count of False, default to 0 if not present

            loi_cca = loi_counts_cca.get(True, 0)  # Get count of True, default to 0 if not present
            nonLoi_cca = loi_counts_cca.get(False, 0)  # Get count of False, default to 0 if not present

            ## DOI
            doi_get = doi_counts.get(True, 0)  # Get count of True, default to 0 if not present
            non_doi = doi_counts.get(False, 0)  # Get count of False, default to 0 if not present


            # Calculate total counts
            total_caa = loi_caa + nonLoi_caa
            total_cca = loi_cca + nonLoi_cca
            total_doi = doi_get + non_doi

            # Calculate percentages
            percentage_caa.append( (loi_caa / total_caa) * 100 if loi_caa > 0 else 0)
            percentage_cca.append((loi_cca / total_cca) * 100 if loi_cca > 0 else 0)
            percentage_doi.append((doi_get / total_doi) * 100 if doi_get > 0 else 0)

        
        df=pd.DataFrame()
        df["Sample"]=campioni
        df["Instability_Index"]=instInd
        df["Expansion_Index"]=expInd
        df["Allele_Ratio"]=histogramRatio
        df["LOI_CAA"]=percentage_caa   
        df["LOI_CCA"]=percentage_cca
        df["CAG_repeatsPeak_Allele_1"]=cag_max_alleles_1
        df["CAG_repeatsPeak_Allele_2"]=cag_max_alleles_2
        df["Max_CAG_observed"]=observed_maxCAG
    
        df.to_excel(output,index=False)
    
        return df

    ################################################## RUN ###################################################################################################################################################################

    # Mostra il logo all'avvio del programma
    print_logo()

    print("Starting strmie")
    print()

    # PIPELINE COMPLETA
    if args.mode == "Complete_Pipeline":

        print("Parameters:")
        print("input: "+input_raw_reads)
        print("resolved input files: "+str(len(input_files)))
        print("output directory: "+path)

        if biological_cutpoint:
            print("cutpoint-based: "+str(biological_cutpoint))
            print("cutpoint: "+str(cutpoint))
        elif cwt:
            print("cwt-based: "+str(cwt))
            print("width: "+str(ampiezza))
            print("interval: "+str(intorno))
        else:
            ## adding nanopore
            if nanopore_mode:
                print("nanopore-mode: " + str(nanopore_mode))
                print("max-roi=" + str(args.np_max_roi))
                print("max-edits=" + str(args.np_max_edits))
                print("max-edits-left=" + str(args.np_max_edits_left))
                print("max-edits-right=" + str(args.np_max_edits_right))
                print("seed-len=" + str(args.np_seed_len))
                print("bestmatch=" + str(args.np_bestmatch))
                print("min-read-len=" + str(args.np_min_read_len))
                print("min-cag-pct=" + str(args.np_min_cag_pct))
                print("cag-pct-scope=" + str(args.np_cag_pct_scope))
                print("allow-caa=" + str(args.np_allow_caa))
            else:
            ## end adding nanopore
    
                print("default-mode")
                print("minimum CAG repeat to consider: "+str(infMin))


        print("CAG-graph: "+str(cag_graph))
        print("CCG-graph: "+str(ccg_andWarning_graph))
        print("threshold Instability Index: "+str(ii_threshold))
        print("threshold Expansion Index: "+str(ei_threshold))
        print()
        print("Start processing...")
        print()

        print("Create dataframe from raw reads")
        #print(input_files)

        list_data=[]

        ## adding nanopore
        #print("Calculate LOI and Freq.")
        #for name in file_names:
        #    path_file=input_raw_reads+name
        #    if c==0:
        #        data=calcola_counts_and_loi(path_file,name)
        #        data["filename"]=name
        #        data=data[data.CAG_repeats>=infMin]
        #        c=c+1
        #    else:
        #        tmp=calcola_counts_and_loi(path_file)
        #        tmp["filename"]=name
        #        tmp=tmp[tmp.CAG_repeats>=infMin]
        #        data=pd.concat([data, tmp])
        
        print("Calculate LOI and Freq.")
        # Each input file is parsed (raw-read regex matching) fully
        # independently of every other file, with no shared state until the
        # concatenation below -- the dominant cost for large sample counts,
        # and a clean multiprocessing target. Parallelized across samples
        # via ProcessPoolExecutor (opt out with -j 1); executor.map()
        # preserves input order regardless of which worker finishes first,
        # so the concatenation order below, and therefore the resulting
        # `data` DataFrame, is identical to the previous purely-sequential
        # version for any given -j.
        nanopore_kwargs = {}
        if nanopore_mode:
            nanopore_kwargs = dict(
                max_roi=args.np_max_roi,
                max_edits=args.np_max_edits,
                max_edits_left=args.np_max_edits_left,
                max_edits_right=args.np_max_edits_right,
                seed_len=args.np_seed_len,
                bestmatch=args.np_bestmatch,
                min_read_len=args.np_min_read_len,
                min_cag_pct=args.np_min_cag_pct,
                cag_pct_scope=args.np_cag_pct_scope,
                allow_caa=args.np_allow_caa,
            )

        tasks = [(path_file, os.path.basename(path_file), nanopore_mode, nanopore_kwargs) for path_file in input_files]
        n_jobs = args.jobs if args.jobs else min(len(tasks), os.cpu_count() or 1, 4)
        n_jobs = max(1, min(n_jobs, len(tasks)))
        print(f"Parallel jobs: {n_jobs}")

        if n_jobs <= 1:
            results = [_process_one_input_file(t) for t in tasks]
        else:
            with ProcessPoolExecutor(max_workers=n_jobs) as executor:
                results = list(executor.map(_process_one_input_file, tasks))

        for i, (name, df_result) in enumerate(results):
            df_result["filename"] = name
            df_result = df_result[df_result.CAG_repeats >= infMin]
            if i == 0:
                data = df_result
            else:
                data = pd.concat([data, df_result])

        ## end adding nanopore
        
        print("Filtering all CAG repetitions lower than 7")
        data=data[data.CAG_repeats>=infMin]



        campioni=list(data.filename.unique())

        dir1="CAG_graphs"
        dir2="CCG_alleles_graphs"
        dir3="warning_case"
        dir4="forced_graphs"
        dir5="raw_counts"


        create1 = os.path.join(path, dir1)
        create2 = os.path.join(path, dir2)
        create3 = os.path.join(path, dir3)
        create4 = os.path.join(path, dir4)
        create5 = os.path.join(path, dir5)

        #createFolder=True

        if cag_graph & ccg_andWarning_graph:
            folders=[create1,create2,create3,create4,create5]
        elif cag_graph:
            folders=[create1,create5]
        elif ccg_andWarning_graph:
            folders=[create2,create3,create4,create5]
        else:
            folders=[create5]

        #    createFolder=False

        #if createFolder:
        for c in folders:
            try:
                os.mkdir(c)
                print("Directory '%s' created" % c)
            except FileExistsError:
                print("Directory '%s' already exists" % c)


        create1=create1+"/"
        create2=create2+"/"
        create3=create3+"/"
        create4=create4+"/"
        create5=create5+"/"


        if cag_graph:
            print("plotting the cag graphs")
            barplot_alleli_samples(data,create1)

        # raw counts
        print("Writing raw counts files")
        for s in list(data.filename.unique()):
            tmp_counts=data[data.filename==s]
            ### salvo dataframe per fare l'istogramma con html report
            tmp_counts.to_csv(create5+str(s)+".csv",index=False)



        out1="report_0.xlsx"
        print("Calculate cag-ccg content, indices and make draft report")
        pear_dataframe=report_to_excel(data,campioni,path+out1,cwt)

        print("Calculate ccg content")
        final,reRun=ccg_count(data,pear_dataframe,create2,create3,ccg_andWarning_graph)
        #print(final.keys())

        if reRun.empty:
            outFile="Final_report.xlsx"
            print("save final report")
            final.to_excel(path+outFile,index=False)
            print("Done")
        else:
            campioni_reRun=list(reRun.filename.unique())

            file_forcedReRun="forced search.xlsx"
            print("force searching of WARNING peaks")
            reRun_forced=force_findingPeaks(reRun,campioni_reRun,path+file_forcedReRun)

            print("calculate ccg content of WARNING peaks")
            final_reRun2,reRun_2=ccg_count(reRun,reRun_forced,create4,create4,ccg_andWarning_graph)
            
            print("preparing final report")
            integrazione=final[final.CAG_repeatsPeak_Allele_1!="warning"]

            integrazione=pd.concat([integrazione,final_reRun2])

            final.index=final.Sample
            integrazione.index=integrazione.Sample

            final.loc[integrazione.index, :] = integrazione[:]

            outFile="Final_report.xlsx"
            print("Save final report")
            final.to_excel(path+outFile,index=False)
    
        ## genera html file
        create_html(path,final,data,cutpoint=cutpoint)
        write_histogram_spreadsheet(path,data)



    # ONLY INDEX CALCULATION
    elif args.mode == "Index_Calculation":

        print("cutpoint: "+str(cutpoint))

        cag_file=pd.read_excel(index_path)

####################################### NEW part 
        # 2) Read and concat CSV in path/raw_counts
        import glob

        raw_counts_dir = os.path.join(path, "raw_counts")
        if not os.path.isdir(raw_counts_dir):
            raise FileNotFoundError(f"Cartella 'raw_counts' non trovata in: {raw_counts_dir}")

        csv_files = sorted(glob.glob(os.path.join(raw_counts_dir, "*.csv")))
        if not csv_files:
            raise FileNotFoundError(f"Nessun .csv trovato in {raw_counts_dir}")

        frames = []
        for fcsv in csv_files:
            dfc = pd.read_csv(fcsv)
            # se manca 'filename' nel csv, derivala dal nome file .csv
            if "filename" not in dfc.columns:
                base = os.path.basename(fcsv)
                if base.endswith(".csv"):
                    base = base[:-4]
                dfc["filename"] = base
            frames.append(dfc)

        # questo è il "data" da usare da qui in avanti
        data = pd.concat(frames, ignore_index=True)

#########################################################################################
    
        df_merged = pd.merge(cag_file, data, left_on="Sample", right_on="filename", how="left")
        campioni=list(df_merged.filename.unique())
        out_indices="indices_calculation.xlsx"
        print("calculate cag-ccg content, indices and make report")
        calculate_indices_fromFile(df_merged,campioni,path+out_indices,cutpoint=cutpoint,ii_threshold=ii_threshold,ei_threshold=ei_threshold)

        print("Done")
 
    print()
    print("The job is done, Thanks for using strmie")

if __name__ == "__main__":
    main()
