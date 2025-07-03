#!/usr/bin/env python
# coding: utf-8

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import os
import argparse
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

def barplot_alleli(df,titolo,name):

    fig = plt.figure(figsize =(13, 10))
    df['CAG_repeats'].value_counts().sort_index().plot.bar()#plt.ylim(0,0.035)
    #plt.ylim(0,t[name].max()+t[name].max()/10)
    plt.rc('xtick', labelsize=7);plt.rc('ytick', labelsize=10);plt.xticks(rotation=90)
    plt.title(titolo)
    plt.xlabel("Number of CAG repeats")
    plt.ylabel("Count the number of repetitions")
    
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
    plt.xlabel("Number of CCG repeats")
    plt.ylabel("Count the number of repetitions")
    
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
##                "strmie - Analyze Huntington's Disease genetic markers from sequencing data"               ##
##                                                                                                            ##
################################################################################################################
    ''' + Style.RESET_ALL)