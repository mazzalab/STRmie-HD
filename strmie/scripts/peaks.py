#!/usr/bin/env python
# coding: utf-8

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.signal as signal
import math
import os
import argparse
from scipy import signal
from scipy.signal import find_peaks

from strmie.scripts.utility import *


def _two_highest_peaks(counts, intorno, min_prominence_ratio=0.025, min_edge_prominence_ratio=0.3):
    """Return the two highest peaks (CAG values) of a value_counts Series.

    Candidate peaks are ranked by topographic prominence (height above the
    deepest valley separating them from a taller peak), not by raw read
    count, computed over the dense per-CAG-value histogram (gaps filled
    with 0 reads). Two kinds of candidates are considered:

    - Interior local maxima, found with scipy.signal.find_peaks(), which
      correctly handles plateaus/ties (common in low-coverage regions,
      where several adjacent CAG values can share the same read count).
    - The two boundary values of the observed CAG range, checked by hand,
      since scipy.signal.find_peaks() requires lower neighbours on BOTH
      sides and is therefore blind to a true peak sitting at the very edge
      of the data (e.g. a large, well-resolved expansion with no reads
      beyond it). Their prominence is computed the same topographic way,
      using only the one side that has real data.

    The tallest candidate is peak_a. Among the rest, peak_b is the
    highest-prominence candidate that is at least `intorno` CAG units away
    from peak_a (closer local maxima are stutter/noise shoulders of the
    same allele) AND whose prominence reaches a minimum fraction of
    peak_a's prominence; otherwise peak_a is reported twice.

    Interior and boundary candidates need different minimum fractions.
    A boundary candidate's prominence is measured against whatever it
    finds while scanning inward until a taller value turns up (or the far
    end of the data, if none does), so a boundary point sitting above a
    long, gently declining stretch of unrelated low-CAG noise can look
    almost as prominent as a real peak; `min_edge_prominence_ratio` (0.3)
    guards against that. Interior candidates cannot inflate their
    prominence this way (both sides are real, independently measured
    data), so a much lower `min_prominence_ratio` (0.025) is enough to
    tell a true but low-coverage second allele (validated down to ~3.5%
    of the tallest peak's prominence in real Dataset 4 ONT data, where PCR
    bias against longer alleles can leave very little read support) apart
    from noise.
    """
    if counts.empty:
        return None

    idx = counts.index
    lo, hi = int(idx.min()), int(idx.max())
    heights = counts.to_dict()
    arr = np.array([heights.get(v, 0) for v in range(lo, hi + 1)], dtype=float)

    candidates = {}  # cag_value -> (height, prominence, is_edge)

    if len(arr) >= 3:
        peak_pos, props = find_peaks(arr, prominence=0)
        for p, prom in zip(peak_pos, props["prominences"]):
            candidates[lo + int(p)] = (arr[p], prom, False)

    def edge_prominence(is_left):
        h0 = arr[0] if is_left else arr[-1]
        scan = range(1, len(arr)) if is_left else range(len(arr) - 2, -1, -1)
        running_min = h0
        for i in scan:
            if arr[i] > h0:
                return h0 - running_min
            running_min = min(running_min, arr[i])
        return h0 - running_min

    if len(arr) >= 2:
        if arr[0] > arr[1] and arr[0] > 0 and lo not in candidates:
            candidates[lo] = (arr[0], edge_prominence(True), True)
        if arr[-1] > arr[-2] and arr[-1] > 0 and hi not in candidates:
            candidates[hi] = (arr[-1], edge_prominence(False), True)
    elif len(arr) == 1 and arr[0] > 0:
        candidates[lo] = (arr[0], arr[0], True)

    if not candidates:
        peak_a = int(counts.idxmax())
        return peak_a, peak_a

    ranked = sorted(candidates.items(), key=lambda kv: kv[1][1], reverse=True)
    peak_a = ranked[0][0]
    prom_a = ranked[0][1][1]

    def acceptable(v, prominence, is_edge):
        if abs(v - peak_a) < intorno:
            return False
        threshold = min_edge_prominence_ratio if is_edge else min_prominence_ratio
        return prominence >= threshold * prom_a

    peak_b = next(
        (v for v, (_, prominence, is_edge) in ranked[1:] if acceptable(v, prominence, is_edge)),
        peak_a,
    )

    return tuple(sorted([peak_a, peak_b]))


def fine_maxPeak_hist_generated_bycutPoint(df,cutpoint):

    tmp=df['CAG_repeats'].value_counts().sort_index().to_frame()
    tmp=tmp.rename({'count':'height_peak'},axis=1)
    #### (1) filtro tutte le ripetizioni CAG minori di 1
    tmp["CAG_repeats"]=list(tmp.index.values)
    tmp=tmp[tmp["CAG_repeats"]>=7]

    if tmp.empty:
        print("ERROR:")
        print("The coverage of the sample: "+str(df.filename.values[0])+" is not sufficient to perform the analysis.") 
        print("#######")
        raise ValueError("Remove it from the folder and run again strmie: "+str(df.filename.values[0]))

    healthy_allele=tmp[tmp["CAG_repeats"]<=cutpoint]
    phenotype_allele=tmp[tmp["CAG_repeats"]>cutpoint]

    #### vedo se esistono questi valori di CAG come cutpoint
    if (phenotype_allele.empty):
        cag1=healthy_allele[healthy_allele["height_peak"]==healthy_allele["height_peak"].max()]["CAG_repeats"].iloc[0]
        cag2=cag1

    else: # caso ideale
        cag1=healthy_allele[healthy_allele["height_peak"]==healthy_allele["height_peak"].max()]["CAG_repeats"].iloc[0]
        cag2=phenotype_allele[phenotype_allele["height_peak"]==phenotype_allele["height_peak"].max()]["CAG_repeats"].iloc[0]

    return cag1,cag2




def find_peaks_two_alleles(df,ampiezza=[5,6,7,8,9,10], intorno=5):

    counts=df['CAG_repeats'].value_counts()
    counts=counts[counts.index>=7]

    if counts.empty:
        print("ERROR_1")
        print("The coverage of the sample: "+str(df.filename.values[0])+" is not sufficient to perform the analysis.")
        print("#######")
        raise ValueError("Remove it from the folder and run again strmie: "+str(df.filename.values[0]))

    return _two_highest_peaks(counts, intorno)



def force_search(t,intorno=5): # t corrisponde al data_campione presente nella funzione report_to_excel

    counts=t.CAG_repeats.value_counts()
    counts=counts[counts.index>=7] ## cambiato da 3 a 10

    # Controllare se è vuoto
    if counts.empty:
        print("WARNING, Sample:")
        print(t["filename"].unique())
        print("No CAG repeats found")
        return 0,0

    result=_two_highest_peaks(counts, intorno)

    return result



def cag_peaks(df, colonna="CAG_repeats",intorno=5):

    counts = df[colonna].value_counts()
    counts = counts[counts.index>=7]

    result = _two_highest_peaks(counts, intorno)

    if result is None:
        return "warning", "warning"

    return result
