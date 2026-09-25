#!/usr/bin/env python
# coding: utf-8

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import os
import sys
import argparse
import shutil
import subprocess
from colorama import Fore, Style
import gzip

def leggi_fasta_gz(file):
    with gzip.open(file, 'rt') as f:  # 'rt' = read text mode
        lines = f.readlines()
    
    ids = []
    seqs = []
    seq = ""

    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            if seq:  # Salva la sequenza precedente
                seqs.append(seq)
                seq = ""
            ids.append(line)
        else:
            seq += line  # Concatena le sequenze su più righe
    
    if seq:
        seqs.append(seq)

    return pd.DataFrame({'ID': ids, 'Seq': seqs})



def leggi_fastq_gz(file):
    df=pd.DataFrame(pd.read_csv(file, sep='\t', header=None,compression='gzip').values.reshape(-1, 4), columns=['ID', 'Seq', 'Sep', 'Qual'])
    return df

def leggi_nomi_file_inDirectory(dir_path): ## solo file con estenzione fastq.gz
    res = []

    for path in os.listdir(dir_path):
        if path.endswith(".fastq.gz"):
            if os.path.isfile(os.path.join(dir_path, path)):
                res.append(path)
        elif path.endswith(".fasta.gz"):
            if os.path.isfile(os.path.join(dir_path, path)):
                res.append(path)
        else:
            print("skip file: "+str(path))

    if len(res) == 0:
        raise TypeError("The file format should be .fastq.gz")

    return res

def resolve_input_paths(input_paths):
    """Resolve user-supplied -f/--input arguments into a flat list of full paths
    to .fastq.gz/.fasta.gz files. Each entry in input_paths may be a directory
    (all matching files inside are included) or a single file."""
    resolved = []

    for p in input_paths:
        if os.path.isdir(p):
            dir_path = p if p.endswith("/") else p + "/"
            for name in leggi_nomi_file_inDirectory(dir_path):
                resolved.append(os.path.join(p, name))
        elif os.path.isfile(p):
            if p.endswith(".fastq.gz") or p.endswith(".fasta.gz"):
                resolved.append(p)
            else:
                raise TypeError(f"Unsupported file format for input file: {p}. Expected .fastq.gz or .fasta.gz")
        else:
            raise FileNotFoundError(f"Input path not found: {p}")

    if len(resolved) == 0:
        raise TypeError("No .fastq.gz or .fasta.gz files found in the specified input path(s).")

    return resolved


# Recognized paired-end filename marker conventions, checked immediately
# before the .fastq.gz/.fasta.gz extension (e.g. sample_R1.fastq.gz /
# sample_R2.fastq.gz, or sample_1.fastq.gz / sample_2.fastq.gz).
_PAIRED_END_MARKERS = [("_R1", "_R2"), ("_1", "_2")]


def _paired_end_marker(basename):
    """If basename ends with a recognized R1/R2 marker right before its
    extension, return (sample_prefix, mate_suffix, ext, is_r1). Else None."""
    for ext in (".fastq.gz", ".fasta.gz"):
        if not basename.endswith(ext):
            continue
        stem = basename[: -len(ext)]
        for r1_suf, r2_suf in _PAIRED_END_MARKERS:
            if stem.endswith(r1_suf):
                return stem[: -len(r1_suf)], r2_suf, ext, True
            if stem.endswith(r2_suf):
                return stem[: -len(r2_suf)], r1_suf, ext, False
    return None


def detect_paired_end_pairs(file_paths):
    """Detect R1/R2 pairs among file_paths by filename convention (same
    directory, matching prefix, recognized R1/R2 marker). Returns
    (pairs, singles): pairs is a list of (sample_prefix, r1_path, r2_path);
    singles is the list of paths not part of a detected pair, unchanged."""
    by_dir = {}
    for p in file_paths:
        by_dir.setdefault(os.path.dirname(p), {})[os.path.basename(p)] = p

    pairs = []
    consumed = set()
    for p in file_paths:
        if p in consumed:
            continue
        d, b = os.path.dirname(p), os.path.basename(p)
        info = _paired_end_marker(b)
        if info is None:
            continue
        prefix, mate_suffix, ext, is_r1 = info
        mate_path = by_dir.get(d, {}).get(prefix + mate_suffix + ext)
        if mate_path is None or mate_path in consumed or mate_path == p:
            continue
        r1_path, r2_path = (p, mate_path) if is_r1 else (mate_path, p)
        pairs.append((prefix, r1_path, r2_path))
        consumed.add(p)
        consumed.add(mate_path)

    singles = [p for p in file_paths if p not in consumed]
    return pairs, singles


def _parse_pear_stats(pear_stdout):
    """Extract Assembled/Discarded/Not-assembled percentages from PEAR's stdout."""
    out = {"Assembled reads (%)": None, "Discarded reads (%)": None, "Not assembled reads (%)": None}
    patterns = {
        "Assembled reads (%)": r"Assembled reads[ .]*:\s*[\d,]+\s*/\s*[\d,]+\s*\(([\d.]+)%\)",
        "Discarded reads (%)": r"Discarded reads[ .]*:\s*[\d,]+\s*/\s*[\d,]+\s*\(([\d.]+)%\)",
        "Not assembled reads (%)": r"Not assembled reads[ .]*:\s*[\d,]+\s*/\s*[\d,]+\s*\(([\d.]+)%\)",
    }
    for key, pat in patterns.items():
        m = re.search(pat, pear_stdout)
        if m:
            out[key] = float(m.group(1))
    return out


def _resolve_pear():
    """Find the pear executable. Checks PATH, then falls back to the bin/
    directory of the currently running Python's own environment -- keeps
    --merge_paired_end (and --bam/--cram, which relies on it internally for
    extracted paired reads) working even when the strmie entry point was
    invoked directly rather than via an activated conda environment."""
    found = shutil.which("pear")
    if found:
        return found
    candidate = os.path.join(sys.prefix, "bin", "pear")
    if os.path.isfile(candidate):
        return candidate
    return None


def merge_paired_end_with_pear(pairs, workdir, min_overlap=10):
    """Merge each (sample_prefix, r1, r2) in pairs with PEAR. Returns
    (merged_paths, stats): merged_paths is {sample_prefix: gzipped assembled
    fastq path}; stats is a list of per-sample PEAR merge statistics dicts.
    Raises RuntimeError if the 'pear' executable can't be found, or if PEAR
    fails for any pair."""
    pear_bin = _resolve_pear()
    if pear_bin is None:
        names = ", ".join(prefix for prefix, _, _ in pairs)
        raise RuntimeError(
            "--merge_paired_end detected R1/R2 file pairs for sample(s): " + names + ", "
            "but the 'pear' executable was not found. Install it with "
            "`conda install -c bioconda pear`, or merge the pairs yourself "
            "(see the Paired-End case study in the documentation) and pass "
            "the merged file per sample to -f/--input instead."
        )

    merged_dir = os.path.join(workdir, "paired_end_merged")
    os.makedirs(merged_dir, exist_ok=True)

    merged_paths = {}
    stats = []
    for prefix, r1, r2 in pairs:
        out_prefix = os.path.join(merged_dir, prefix)
        cmd = [pear_bin, "-f", r1, "-r", r2, "-v", str(min_overlap), "-o", out_prefix]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(
                f"PEAR failed to merge paired-end reads for sample '{prefix}':\n"
                f"{result.stdout}\n{result.stderr}"
            )

        assembled = out_prefix + ".assembled.fastq"
        if not os.path.isfile(assembled):
            raise RuntimeError(
                f"PEAR did not produce an assembled output for sample '{prefix}' "
                f"(expected {assembled})."
            )

        # Named after the shared sample prefix (not "<prefix>.assembled.fastq.gz")
        # so the merged pair reports under a clean sample name downstream.
        gz_path = os.path.join(merged_dir, prefix + ".fastq.gz")
        with open(assembled, "rb") as f_in, gzip.open(gz_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(assembled)
        merged_paths[prefix] = gz_path

        sample_stats = _parse_pear_stats(result.stdout)
        sample_stats["Sample"] = prefix
        stats.append(sample_stats)

    return merged_paths, stats

def barplot_alleli(df,titolo,name):

    fig = plt.figure(figsize =(13, 10))
    df['CAG_repeats'].value_counts().sort_index().plot.bar()#plt.ylim(0,0.035)
    #plt.ylim(0,t[name].max()+t[name].max()/10)
    plt.rc('xtick', labelsize=7);plt.rc('ytick', labelsize=10);plt.xticks(rotation=90)
    plt.title(titolo)
    plt.xlabel("CAG repeats")
    plt.ylabel("Read Counts")
    
    plt.savefig(name)
    plt.close(fig)


def barplot_alleli_samples(df,path_outIMG):
    samples=list(df.filename.unique())
    for s in samples:
        tmp=df[df.filename==s]
        barplot_alleli(tmp,"CAG content of Alleles",path_outIMG+s+".jpg")
        ### salvo dataframe per fare l'istogramma
        #tmp.to_excel(path_outIMG+s+".xlsx")

    barplot_alleli(df,"CAG content of Alleles for all Samples",path_outIMG+"All_samples"+".jpg")

def barplot_alleli_ccg(df,titolo,name):
    fig = plt.figure(figsize =(13, 10))
    df['CCG_repeats'].value_counts().sort_index().plot.bar()

    plt.rc('xtick', labelsize=7);plt.rc('ytick', labelsize=10);plt.xticks(rotation=90)
    plt.title(titolo)
    plt.xlabel("CCG repeats")
    plt.ylabel("Read Counts")
    
    plt.savefig(name)
    plt.close(fig)

def create_df_distribution(df):
    df_pre_norm=df['CAG_repeats'].value_counts().to_frame()
    df_pre_norm=df_pre_norm.rename({'CAG_repeats':'height_peak'},axis=1)
    df_pre_norm["CAG_repeat"]=df_pre_norm.index
    return df_pre_norm




# Funzione per stampare il logo
def print_logo():
    print(Fore.GREEN + Style.BRIGHT + '''
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
    ''' + Style.RESET_ALL)
