#!/usr/bin/env python
# coding: utf-8

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

from strmie.scripts.utility import *






def instabilityIndex(df_pre_norm,cag1,cag2,pcrFiltering=False):

    ### Only reject the "too close/ambiguous to order" case (cag2 == cag1+1
    ### or cag1 > cag2). cag1 == cag2 is a different, legitimate situation:
    ### peak-calling found only one allele (homozygous call, or a second
    ### allele that didn't clear detection), and STRmie-HD reports that
    ### single peak as both cag1 and cag2 by convention (see
    ### _two_highest_peaks). That single peak is a perfectly valid anchor
    ### for II, not an error condition, so it must not hit this guard.
    if cag1!=cag2 and cag1+1>=cag2: ### Warning
        instability_index="It cannot be calculated"
    else:
        ### prendo solo il secondo picco ##NEW
        if cag1==cag2:
            ### No second allele exists to bound the valley search on the
            ### left, so use the full observed range up to the single
            ### called peak instead of [cag1, cag2] collapsing to one bin
            ### -- otherwise the contraction/left side of the distribution
            ### around the sole peak would be silently dropped entirely.
            tmp=df_pre_norm[df_pre_norm["CAG_repeat"]<=cag2]
        else:
            tmp=df_pre_norm[df_pre_norm["CAG_repeat"]>=cag1]
            tmp=tmp[tmp["CAG_repeat"]<=cag2]
        minima_h=tmp["count"].min()
        global line_h
        ### When multiple CAG bins tie at the minimum count (common in sparse/noisy
        ### regions), take the one closest to cag2 rather than an arbitrary tied
        ### bin: .iloc[0] picked whichever tied row happened to come first in the
        ### unsorted per-count table, which is not deterministic across runs and
        ### could anchor the summation window at a bin far from the true valley.
        line_h = tmp.loc[tmp['count'] == minima_h, 'CAG_repeat'].max()
        df_pre_norm=df_pre_norm[df_pre_norm["CAG_repeat"]>=line_h]
        ### NEW

        ### normalizzo per la somma delle peak height
        df_pre_norm=df_pre_norm.rename({'count':'height_peak'},axis=1)
        df_norm=df_pre_norm[['height_peak']].div(df_pre_norm[['height_peak']].sum(axis=0), axis=1)
        df_norm["CAG_repeat"]=df_pre_norm["CAG_repeat"]

        #if (df_norm[df_norm.CAG_repeat == cag1].empty) | (df_norm[df_norm.CAG_repeat == cag2].empty):
            #print("Exception")
        if (df_norm[df_norm.CAG_repeat == cag2].empty):
            return "Necessary data points for CAG1 or CAG2 are missing"

        else:
            if pcrFiltering!=False: ### Filtro PCR 
                #if not df_norm[df_norm.CAG_repeat == cag1].empty and not df_norm[df_norm.CAG_repeat == cag2].empty:
                #maxPeakHeight_allele1 = df_norm['height_peak'][df_norm.CAG_repeat == cag1].values[0] # NEW
                #maxPeakHeight_allele2 = df_norm['height_peak'][df_norm.CAG_repeat == cag2].values[0]
                #else:
                #    return "Necessary data points for CAG1 or CAG2 are missing"
                #if maxPeakHeight_allele1>=maxPeakHeight_allele2:
                 #   less20percent=maxPeakHeight_allele2*0.20
                #else:
                #    less20percent=maxPeakHeight_allele1*0.20
                maxPeakHeight_allele2 = df_norm['height_peak'][df_norm.CAG_repeat == cag2].values[0]
                less20percent=maxPeakHeight_allele2*pcrFiltering # NEW
                df_norm=df_norm[df_norm['height_peak']>less20percent]

            ### i picchi normalizzati vengono moltiplicati per il valore di ordinamento dato dal picco del secondo allele (RANGO)
            pre=len(df_norm[df_norm["CAG_repeat"]<cag2])
            post=len(df_norm[df_norm["CAG_repeat"]>cag2])
            pre_ord= [*range(-pre, 0, 1)]
            post_ord= [*range(0, post+1, 1)]
            ordinamento=pre_ord+post_ord
            df_norm=df_norm.sort_values('CAG_repeat')
            df_norm["ordinamento"]=ordinamento
            df_norm["heightXchanges"]=df_norm['height_peak']*df_norm["ordinamento"]

            ### sommo i volori ottenuti al punto precedente per ottenere l'II
            instability_index=df_norm["heightXchanges"].sum()

    return instability_index

def expansionIndex(df_pre_norm,cag1,cag2,pcrFiltering=False):

    ### Same relaxation as instabilityIndex: cag1 == cag2 means only one
    ### allele was called (homozygous, or no second peak cleared
    ### detection), and that single peak should still be usable as the
    ### reference for EI -- the formula below never actually uses cag1
    ### except for this guard, so cag1==cag2 needs no further special
    ### casing past this point. Only "too close to order" (cag2==cag1+1,
    ### or cag1>cag2) remains an error.
    if cag1!=cag2 and cag1+1>=cag2: ### warning
        expansion_index="It cannot be calculated"
    else:
        ### normalizzo per l'altezza del picco del secondo allele 
        df_pre_norm=df_pre_norm.rename({'count':'height_peak'},axis=1)
        if not df_pre_norm[df_pre_norm.CAG_repeat == cag2].empty:
            maxPeakHeight_allele2 = df_pre_norm['height_peak'][df_pre_norm.CAG_repeat == cag2].values[0]
        else:
            return "Necessary data points for CAG1 or CAG2 are missing"
        df_norm=df_pre_norm[['height_peak']].div(maxPeakHeight_allele2)
        df_norm["CAG_repeat"]=df_pre_norm["CAG_repeat"]

        if pcrFiltering!=False: ### Filtro PCR
            
            norm_maxPeakHeight_allele2 = df_norm['height_peak'][df_norm.CAG_repeat == cag2].values[0]
            less3percent=norm_maxPeakHeight_allele2*pcrFiltering   
            df_norm=df_norm[df_norm["height_peak"]>less3percent]


        ### Prendo soltanto i picchi che si trovano a destra del picco del secondo allele (picco del secondo allele compreso)
        df_norm=df_norm[df_norm['CAG_repeat']>=cag2]

        ### i picchi normalizzati vengono moltiplicati per il valore di ordinamento dato dal picco del secondo allele (RANGO)
        ordinamento=[*range(0, len(df_norm), 1)]
        df_norm=df_norm.sort_values('CAG_repeat')
        df_norm["ordinamento"]=ordinamento
        df_norm["heightXchanges"]=df_norm['height_peak']*df_norm["ordinamento"]
        
        ### sommo i volori ottenuti al punto precedente per ottenere l'EI
        expansion_index=df_norm["heightXchanges"].sum()
    
    return expansion_index



def histogramRatioIndex(df_pre_norm,cutpoint=39):

        ### normalizzo per la somma delle peak height
    df_pre_norm=df_pre_norm.rename({'count':'height_peak'},axis=1)
    df_norm=df_pre_norm[['height_peak']].div(df_pre_norm[['height_peak']].sum(axis=0), axis=1)
    df_norm["CAG_repeat"]=df_pre_norm["CAG_repeat"]

    healthy_allele=df_norm[df_norm["CAG_repeat"]<=cutpoint]
    phenotype_allele=df_norm[df_norm["CAG_repeat"]>cutpoint]

    healthy_area=healthy_allele["height_peak"].sum()
    phenotype_area=phenotype_allele["height_peak"].sum()
    ratio_index=phenotype_area/healthy_area

    return round(ratio_index,2)
