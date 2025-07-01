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

import torch
from torch import nn, optim, distributions
from torch.utils.data import Dataset, DataLoader, ConcatDataset
from torchvision import transforms
import pickle as pk
from scipy import stats
from sklearn.model_selection import train_test_split

from strmie.scripts.utility import *
from strmie.scripts.html_generator import *
from strmie.scripts.custom_transformers import *
from strmie.scripts.fully_conv_NN import FindPeaksModel
from strmie.scripts.indices import *
from strmie.scripts.pattern import *
from strmie.scripts.peaks import *

################################################### CNN #################################################################################
def coppie_PE(lista_files):
    lista_files.sort()
    forward=[]
    reverse=[]
    for i,f in enumerate(lista_files):
        if (i+1) % 2 == 1:
            forward.append(f)

        if (i+1) % 2 == 0:
            reverse.append(f)
    
    print("Check if they are PE and if they are all coupled...")

    if len(forward)==len(reverse):
        for i,e in enumerate(forward):
            diff_letters = sum( e[j] != reverse[i][j] for j in range(len(e)) )
            if diff_letters>1:
                print("WARNING: some files do not have the same name for reverse and forward")
    else:
        print("WARNING: some files are not paired")
        raise ValueError("Some files are not paired")
    print("Ok, they are PE")

    return forward, reverse


class UnknownSeriesDataset(Dataset):

  def __init__(self, all_ids, transform=None, max_len=None):
    self.all_ids = all_ids
    self.transform = transform
    self.max_len = max_len

  def __len__(self):
    return len(self.all_ids)

  def __getitem__(self, i):
    X = torch.tensor(hist_data.loc[self.all_ids[i]].values, dtype=torch.float32)
    y = torch.tensor([torch.nan, torch.nan], dtype=torch.float32)
    if not self.max_len is None:
      X = X[:,:self.max_len]
    if not self.transform is None:
        X, y = self.transform([X, y])
    return X.type(torch.float32), y.type(torch.float32)


def cnn_find_peaks(hist_data):
    #print(hist_data)
    INPUT_LEN = 81
    input_ids = hist_data.index.get_level_values(0)[::2]

    trivial_transformer = transforms.Compose([Log1p(),ClipLabels(lower=1, upper=INPUT_LEN),MinMaxScaling(lower=0, upper=1),])

    #print(input_ids)
    input_dataset = UnknownSeriesDataset(input_ids, transform=trivial_transformer, max_len=INPUT_LEN)
    input_dataloader = DataLoader(input_dataset, batch_size=1, shuffle=False, num_workers=2)

    model = pk.load(open("/data3/analysis3/strmie/Scripts/model_object.pk", "rb"))
    model.load_state_dict(torch.load("state_dict.weight", map_location="cpu"))

    predicted_peaks = pd.DataFrame([], columns=["peak_1", "peak_2"])
    #peaks_probabilities = pd.DataFrame(
        #[], columns=pd.MultiIndex.from_product([["peak_1", "peak_2"], np.arange(INPUT_LEN)])
    #)
    for i, (X, _) in enumerate(input_dataloader):
        predicted_peaks.loc[input_ids[i], ["peak_1", "peak_2"]] = model.find_peaks(X)[0].detach().numpy()
        #peaks_probabilities.loc[input_ids[i], :] = model(X)[0].detach().numpy().ravel()
        #tmp=model(X)[0].detach().numpy()
        #peaks_probabilities.loc[input_ids[i], "peak_1"]=tmp[0]
        #peaks_probabilities.loc[input_ids[i], "peak_2"]=tmp[1]

    return predicted_peaks



def metrics_for_predicted_peaks(data,campioni,output,predicted_peaks):

    instInd=[]
    expInd=[]
    cag_max_alleles_1=[]
    cag_max_alleles_2=[]
    percentage=[]
    histogramRatio=[]

    data["sample"]=data.filename.apply(lambda x: x.split("_L001")[0])
    predicted_peaks["sample"]=list(predicted_peaks.index)
    for c in campioni:
        data_campione=data[data.filename==c]
        sample=data["sample"][data.filename==c].values[0]
        cag_max_1=predicted_peaks["peak_1"][predicted_peaks["sample"]==sample].values[0]
        cag_max_2=predicted_peaks["peak_2"][predicted_peaks["sample"]==sample].values[0]

        df_distrib=create_df_distribution(data_campione)
        ii=instabilityIndex(df_distrib,cag_max_1,cag_max_2)  # instability Index , ti
        ei=expansionIndex(df_distrib,cag_max_1,cag_max_2)    # expansion index ,te
        instInd.append(ii)
        expInd.append(ei)
        cag_max_alleles_1.append(cag_max_1)
        cag_max_alleles_2.append(cag_max_2)
        histogramRatio.append(histogramRatioIndex(df_distrib),cutpoint)

        # Data filtering
        df_loi = data_campione[data_campione.CAG_repeats >= cag_max_2]

            # Count LOI values
        loi_counts_caa = df_loi.LOI_CAA.value_counts()
        loi_counts_cca = df_loi.LOI_CCA.value_counts()

            # Check for True and False in the index
        loi_caa = loi_counts_caa.get(True, 0)  # Get count of True, default to 0 if not present
        nonLoi_caa = loi_counts_caa.get(False, 0)  # Get count of False, default to 0 if not present

        loi_cca = loi_counts_cca.get(True, 0)  # Get count of True, default to 0 if not present
        nonLoi_cca = loi_counts_cca.get(False, 0)  # Get count of False, default to 0 if not present

            # Calculate percentages based on the presence of LOI values
        if loi_caa > 0 and nonLoi_caa > 0:  # Some have LOI
            percentage_caa.append((loi_caa / (loi_caa + nonLoi_caa)) * 100)
        elif loi_caa > 0 and nonLoi_caa == 0:  # All have LOI
            percentage_caa.append(100)
        else:  # None have LOI
            percentage_caa.append(0)

            # Calculate percentages based on the presence of LOI values
        if loi_cca > 0 and nonLoi_cca > 0:  # Some have LOI
            percentage_cca.append((loi_cca / (loi_caa + nonLoi_cca)) * 100)
        elif loi_cca > 0 and nonLoi_cca == 0:  # All have LOI
            percentage_cca.append(100)
        else:  # None have LOI
            percentage_cca.append(0)

        
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
    
    df.to_excel(str(output))
    
    return df
