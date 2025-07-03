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




def find_peaks_two_alleles(df,ampiezza=[5,6,7,8,9,10], intorno=6):

    tmp=df['CAG_repeats'].value_counts().sort_index().to_frame()
    tmp=tmp.rename({'count':'height_peak'},axis=1)
    #### (1) filtro tutte le ripetizioni CAG minori di 1
    tmp["CAG_repeats"]=list(tmp.index.values)
    tmp=tmp[tmp["CAG_repeats"]>=7]

    if tmp.empty:
        print("ERROR_1")
        print("The coverage of the sample: "+str(df.filename.values[0])+" is not sufficient to perform the analysis.") 
        print("#######")
        raise ValueError("Remove it from the folder and run again strmie: "+str(df.filename.values[0]))

    altezze=tmp["height_peak"].values
    peak_indices= signal.find_peaks_cwt(altezze,widths=ampiezza)
    
    if len(peak_indices)==2: 
        peaks=[int(tmp[tmp.height_peak==altezze[peak_indices[0]]].index[0]),int(tmp[tmp.height_peak==altezze[peak_indices[1]]].index[0])] ## ripetizioni CAG (non altezze)
        df_check_peak_1=tmp.loc[(tmp.CAG_repeats<(peaks[0]+intorno)) & (tmp.CAG_repeats>(peaks[0]-intorno))]
        check_cag_1=df_check_peak_1.CAG_repeats[df_check_peak_1.height_peak==df_check_peak_1.height_peak.max()].values[0]

        if (peaks[0]<40) & (peaks[1]>26):
            df_check_peak_2=tmp.loc[(tmp.CAG_repeats<(peaks[1]+intorno)) & (tmp.CAG_repeats>(peaks[1]-intorno))]
            check_cag_2=df_check_peak_2.CAG_repeats[df_check_peak_2.height_peak==df_check_peak_2.height_peak.max()].values[0]
            return check_cag_1,check_cag_2
        else:
            #print("warning peak")
            return "warning: "+str(peaks[0]),"warning: "+str(peaks[1])
    else:
        return "warning","warning"



def force_search(t,intorno=6): # t corrisponde al data_campione presente nella funzione report_to_excel
    name=t.filename
    tmp_name=t["filename"].unique()
    
    t=t.CAG_repeats.value_counts().to_frame()
    t["CAG_repeat"]=t.index
    t.columns.names = [None]
    t=t.rename({'count':'height_peak'},axis=1)
    t=t[t.CAG_repeat>=7] ## cambiato da 3 a 10 
    t.sort_values(by=["CAG_repeat"],inplace=True)
    media=t.height_peak.mean()
    massimo=t.height_peak.max()

    # Controllare se è vuoto
    if t.empty:
        print("WARNING, Sample:")
        print(t["filename"].unique())
        print("No CAG repeats found")
        return 0,0


    filtro = t[t.height_peak == massimo]

    # Se ci sono più di una riga nel DataFrame filtrato, prendi l'ultima riga
    if len(filtro) > 1:
        cag_max=filtro.tail(1)['CAG_repeat'].values[0]
    else:
        cag_max=filtro['CAG_repeat'].values[0]  # Se c'è solo una riga, prendi semplicemente quella


    #cag_max=t.CAG_repeat[t.height_peak==massimo].values[0]
    altezze=list(t["height_peak"])

    peaks, properties = find_peaks(altezze, distance=5,wlen=3)#### PARAMETRIZZABILE:The required minimum number of data points between peaks.
    # If there are more than two peaks, select the two with the highest peak heights
    if len(peaks) > 2:
        max_heights=[]

        for p in peaks:
            ## controllo il primo picco
            df_check_peak=t.loc[(t.height_peak<(altezze[p]+intorno)) & (t.height_peak>(altezze[p]-intorno))]
            check_cag=df_check_peak.height_peak.max()

            filtro2=df_check_peak[df_check_peak.height_peak==check_cag]
            # Se ci sono più di una riga nel DataFrame filtrato, prendi l'ultima riga
            if len(filtro2) > 1:
                cag_prova=list(filtro2.tail(1)["height_peak"].values)
            else:
                cag_prova=list(filtro2["height_peak"].values)  # Se c'è solo una riga, prendi semplicemente quella

            #cag_prova=list(df_check_peak["height_peak"][df_check_peak.height_peak==check_cag].values)
            if len(cag_prova)==1:
                max_heights.append(check_cag)
        
        if len(max_heights)==1:
            h_cag1_tmp=max_heights[0]
            h_cag2_tmp=max_heights[0]
        else:
            h_cag1_tmp=max(max_heights)
            max_heights.remove(h_cag1_tmp)
            h_cag2_tmp=max(max_heights)

        cag1_tmp=t["CAG_repeat"][t.height_peak==h_cag1_tmp].values[0]
        cag2_tmp=t["CAG_repeat"][t.height_peak==h_cag2_tmp].values[0]

        if cag1_tmp>cag2_tmp:
            cag1=cag2_tmp
            cag2=cag1_tmp
        else:
            cag1=cag1_tmp
            cag2=cag2_tmp

    elif len(peaks)==1:
        cag1=peaks[0]
        cag2=peaks[0]
    
    elif len(peaks) < 1:
        cag1=0
        cag2=0

    else:
        cag1=t.CAG_repeat[t.height_peak==altezze[peaks[0]]].values[0]
        cag2=t.CAG_repeat[t.height_peak==altezze[peaks[1]]].values[0]


    return cag1,cag2



def cag_peaks(df, colonna="CAG_repeats",intorno=5):

    bins=(int(df[colonna].max()) - int(df[colonna].min()))

    if bins==0:
        return int(df[colonna].max()),int(df[colonna].max())

    # Compute histogram: `hist_values` stores frequency counts, `bin_edges` contains bin edges
    hist_values, bin_edges = np.histogram(df[colonna], bins=bins)

    # Find peaks in the histogram
    # `height=0` ensures only actual peaks (above zero) are considered
    # `distance=5` ensures peaks are separated by at least 5 bins to avoid local maxima within the same peak
    peaks, properties = find_peaks(hist_values, height=0, distance=intorno)
    
    # If fewer than two peaks are detected, issue a warning and exit
    if len(peaks) == 0:
        cag1 = "warning"
        cag2 = "warning"
    elif len(peaks) == 1:
        # Se c'è solo un picco, entrambi i valori sono uguali a quel picco
        single_peak_index = peaks[0]
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2  # Calcolo dei centri dei bin
        cag1 = int(bin_centers[single_peak_index])
        cag2 = cag1
    else:
        # Extract peak heights from detected peaks
        peak_heights = properties["peak_heights"]

        # Sort peaks by their heights and get the indices of the two tallest peaks
        sorted_indices = np.argsort(peak_heights)[-2:]  # Selects last two indices (highest peaks)
        top_two_peaks = peaks[sorted_indices]  # Get corresponding peak positions

        # Compute bin centers from edges for accurate peak positioning
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        # Retrieve and print the two highest peak CAG values
        peak_cag_values = bin_centers[top_two_peaks]

        vettore=[int(peak_cag_values[1]),int(peak_cag_values[0])]
        vettore_ordinato = sorted(vettore)  # Restituisce una nuova lista ordinata
        cag1 = vettore_ordinato[0]
        cag2 = vettore_ordinato[1]
        
    return cag1,cag2
