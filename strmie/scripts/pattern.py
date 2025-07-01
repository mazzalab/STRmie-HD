#!/usr/bin/env python
# coding: utf-8

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import os
import sys

from strmie.scripts.utility import *
from strmie.scripts.indices import *
from strmie.scripts.peaks import *



def htt_exact_match(seq):

    pattern = r"(?P<CAG>(cag)+)((?P<LOI_CAA>caacag)?){0,2}ccg(?P<LOI_CCA>cca)?(?P<CCG>(ccg)+)*(cct)*"

    # Esegue la ricerca nella sequenza 'seq' usando il pattern definito
    # re.MULTILINE permette di trattare la sequenza come se fosse su più righe
    # re.IGNORECASE ignora la differenza tra maiuscole e minuscole
    match = re.search(pattern, seq, re.MULTILINE | re.IGNORECASE)

        # Se non viene trovato nessun match, ritorna 4 valori None
    if match==None:
        #print("did not find")
        return(None,None,None,None,None)
    else:
            # Se viene trovato un match, ottieni i gruppi di cattura
        groups = match.groups()  # Raccoglie tutti i gruppi catturati dal match
        groups_dict = match.groupdict()  # Ottiene un dizionario con i gruppi per nome
        
            # Calcola la lunghezza della sequenza "CAG" (in termini di numero di triplette)
        cag_len = len(match.group("CAG"))/3
        
            # Se il gruppo "CCG" non è presente, imposta la sua lunghezza a 0
        if match.group("CCG")==None:
            ccg_len = 0
        else:
                # Altrimenti calcola la lunghezza di "CCG" in termini di triplette
            ccg_len = len(match.group("CCG"))/3

        # Verifica se il gruppo "LOI" è presente. Se non è presente, significa che non è stato trovato
        # il pattern "caa", quindi 'is_loi' è False. Altrimenti, è True se non è stato trovato il pattern CAA
        is_loi_caa = match.group("LOI_CAA") is None

        # Verifica la presenza di "CCA" (se "cca" è presente in "ccgcca")
        is_loi_cca = match.group("LOI_CCA") is None
        
        # Conta il numero di occorrenze di "caa" nella sequenza identificata
        matched_seq = match.group(0)
        caa_count = matched_seq.lower().count("caa")
        is_doi = True if caa_count == 2 else False

         # Debugging opzionale
        #if is_doi:
            #print("DOI trovato:", is_doi, "\nSequenza:", seq, "\n")
        #if is_loi_caa:
            #print("LOI CAA trovato:", is_loi_caa, "\nSequenza:", seq, "\n")
        #if is_loi_cca and not is_loi_caa:
         #   print("LOI CCA trovato:", is_loi_cca, "\nSequenza:", seq, "\n")
        
        # Restituisce la lunghezza della sequenza "CAG", lo stato del "LOI" e la lunghezza della sequenza "CCG"
        return(cag_len,is_loi_caa,is_loi_cca,ccg_len,is_doi)




def calcola_counts_and_loi(path_file,name=None):
    
    if path_file.endswith('.fastq.gz'):
        df=leggi_fastq_gz(path_file)
    elif path_file.endswith('.fasta.gz'):
        df=leggi_fasta_gz(path_file)
    else:
        print("ERROR: The file format should be .fastq.gz")

    #df=leggi_fastq_gz(path_file)
    #print(df) ## questo df contiene id seq e qualità di tutte le reads per un singolo campione
    df_count=pd.DataFrame()
    
    cag_len=[]
    ccg_len=[]
    loi_caa=[]
    loi_cca=[]
    id_sequenze=[]
    sequenze=[]
    doi=[]

    for seq,idseq in zip(df.Seq, df.ID):    
        cag,is_loi_caa,is_loi_cca,ccg,has_doi = htt_exact_match(seq)
        cag_len.append(cag)
        ccg_len.append(ccg)
        loi_caa.append(is_loi_caa)
        loi_cca.append(is_loi_cca)
        id_sequenze.append(idseq)
        sequenze.append(seq)
        doi.append(has_doi)

    df_count["ID"]=id_sequenze
    df_count["CAG_repeats"]=cag_len
    df_count["CCG_repeats"]=ccg_len
    df_count["LOI_CAA"]=loi_caa
    df_count["LOI_CCA"]=loi_cca
    df_count["Seq"]=sequenze
    df_count["DOI"]=doi
    
    df_count=df_count.loc[((df_count.LOI_CAA==True) | (df_count.LOI_CAA==False))] ## Serve a togliere i warning, non serve farlo anche per cca
    #print(df_count)##questo df è relativo ad un singolo campione con tutte le ripetizioni presenti in ogni reads del campione
    #df_count=df_count[df_count["CAG_repeats"]>1]
    if df_count.empty:
        sys.exit("ERROR: No CAG repeats found for this sample. Please rerun the analysis after removing this sample:"+str(name))

    return df_count





def ccg_count(df,df_report,path_outIMG,path_warn,save):
        
    samples=list(df_report.Sample.unique())
    ccg_counts_allele_1=[]
    ccg_counts_allele_2=[]
    reRun=pd.DataFrame()

    for s in samples: 
        whatIsIt=list(df_report.CAG_repeatsPeak_Allele_1[df_report.Sample==s])
        whatIsIt=whatIsIt[0]
        
        if type(whatIsIt) == str:
            ccg_counts_allele_1.append("warning")
            ccg_counts_allele_2.append("warning")
            if save:
                barplot_alleli(df[df.filename==s],"WARNING: CAG content of sample -noSingleAllele-",path_warn+"WARNING_CAG_"+s+"_sample.jpg")
            reRun=pd.concat([reRun, df[df.filename==s]])
        else:
        ## ALLELE 1
            max_repeatitionCAG=df_report["CAG_repeatsPeak_Allele_1"][df_report.Sample==s].values[0]
            tmp=df.loc[((df.filename==s)&(df.CAG_repeats==max_repeatitionCAG))]
            
            if tmp.empty:
                ccg_counts_allele_1.append(0)
                ccg_counts_allele_2.append(0)
            
            else:
                most_frequentCCG=tmp.CCG_repeats.value_counts().index.values[0]
                if save:
                    barplot_alleli_ccg(tmp,"CCG content Allele 1",path_outIMG+s+"allele_1.jpg")
                ccg_counts_allele_1.append(most_frequentCCG)
            ## ALLELE 2
                max_repeatitionCAG=df_report["CAG_repeatsPeak_Allele_2"][df_report.Sample==s].values[0]
                tmp=df.loc[((df.filename==s)&(df.CAG_repeats==max_repeatitionCAG))]
                most_frequentCCG=tmp.CCG_repeats.value_counts().index.values[0]
                if save:
                    barplot_alleli_ccg(tmp,"CCG content Allele 2",path_outIMG+s+" allele_2.jpg")
                ccg_counts_allele_2.append(most_frequentCCG)

    df_report["CCG_allele_1"]=ccg_counts_allele_1
    df_report["CCG_allele_2"]=ccg_counts_allele_2
    
    return df_report,reRun






def calculate_indices_fromFile(data,campioni,output,cutpoint=27):
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
        cag_max_1=data_campione["CAG_Allele_1"].values[0]
        cag_max_2=data_campione["CAG_Allele_2"].values[0]

        cag_max_alleles_1.append(cag_max_1)
        cag_max_alleles_2.append(cag_max_2)
        

        df_distrib=create_df_distribution(data_campione)
        observed_maxCAG.append(df_distrib["CAG_repeat"].max())
        ii=instabilityIndex(df_distrib,cag_max_1,cag_max_2)  # instability Index , ti
        ei=expansionIndex(df_distrib,cag_max_1,cag_max_2)    # expansion index ,te
        instInd.append(ii)
        expInd.append(ei)
        histogramRatio.append(histogramRatioIndex(df_distrib,cutpoint))
    
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
        percentage_caa.append((loi_caa / total_caa) * 100 if loi_caa > 0 else 0)
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

    df.to_excel(str(output))
    
    return df