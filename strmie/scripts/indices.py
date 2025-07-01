#!/usr/bin/env python
# coding: utf-8

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

from strmie.scripts.utility import *






def instabilityIndex(df_pre_norm,cag1,cag2,pcrFiltering=False):

    if cag1+1>=cag2: ### Warning
        instability_index="It cannot be calculated"
    else:
        ### prendo solo il secondo picco ##NEW
        tmp=df_pre_norm[df_pre_norm["CAG_repeat"]>=cag1]
        tmp=tmp[tmp["CAG_repeat"]<=cag2]
        minima_h=tmp["count"].min()
        global line_h
        line_h = tmp.loc[tmp['count'] == minima_h, 'CAG_repeat'].iloc[0]
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

    if cag1+1>=cag2: ### warning
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
        maxPeakHeight_norm_allele2 = df_norm['height_peak'][df_norm.CAG_repeat == cag2].values[0]
        indice_maxPeakHeight=df_norm["CAG_repeat"][df_norm["height_peak"]==maxPeakHeight_norm_allele2].values[0]
        df_norm=df_norm[df_norm['CAG_repeat']>=indice_maxPeakHeight]

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
