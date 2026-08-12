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


def _two_highest_peaks(counts, intorno, min_prominence=0.2):
    """Return the two highest peaks (CAG values) of a value_counts Series.

    A CAG value is a local maximum if its read count is strictly higher than
    both of its immediate integer neighbours (missing neighbours count as 0).
    This intentionally does NOT use scipy.signal.find_peaks(), which requires
    a peak to have lower neighbours on BOTH sides and is therefore blind to a
    true peak sitting at the edge of the observed CAG range (e.g. a large,
    well-resolved expansion with no reads beyond it).

    The two tallest local maxima are reported as the two alleles, provided
    they are more than `intorno` CAG units apart; local maxima within
    `intorno` of the tallest one are treated as stutter/noise shoulders of
    the same allele and skipped. This keeps genuinely close-but-distinct
    alleles (a few CAG units apart) from being merged into one, while still
    consolidating same-allele sequencing noise into a single call.

    A candidate second peak is only accepted if its height is at least
    `min_prominence` of the tallest peak's height. Without this, the CAG>=7
    floor applied upstream creates an artificial edge at CAG=7 that the
    "missing neighbour counts as 0" rule treats exactly like a genuine
    edge peak, even when it is really just the tail of a low-CAG noise
    population that got truncated by the floor rather than a real allele.
    Real second alleles observed in this project's validated data sit at
    >=49% of the tallest peak's height, while floor-noise artifacts sit
    below 13%, so 20% cleanly separates the two without needing to special
    case the floor value itself.
    """
    if counts.empty:
        return None

    heights = counts.to_dict()

    local_maxima = [
        (int(v), h) for v, h in heights.items()
        if h > heights.get(v - 1, 0) and h > heights.get(v + 1, 0)
    ]

    if not local_maxima:
        peak_a = int(counts.idxmax())
        return peak_a, peak_a

    local_maxima.sort(key=lambda x: x[1], reverse=True)

    peak_a, height_a = local_maxima[0]
    peak_b = next(
        (v for v, h in local_maxima[1:]
         if abs(v - peak_a) >= intorno and h >= min_prominence * height_a),
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
