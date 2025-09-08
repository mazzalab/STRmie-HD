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

    
from scipy import stats

from strmie.scripts.utility import *
from strmie.scripts.html_generator import *
from strmie.scripts.indices import *
from strmie.scripts.pattern import *
from strmie.scripts.peaks import *


from colorama import Fore, Style

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
    ##                "strmie - Analyze Huntington's Disease genetic markers from sequencing data"               ##
    ##                                                                                                            ##
    ################################################################################################################
    ''' + Style.RESET_ALL + Fore.MAGENTA + "usage: strmie.py -f /path/to/raw_reads_file.fastq.gz -o /path/directory_output/ command [subcommand] parameters" + Style.RESET_ALL , add_help=True, formatter_class=argparse.RawDescriptionHelpFormatter)

    # Definizione della modalità
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--mode", choices=["Complete_Pipeline", "Index_Calculation"], help=Fore.CYAN + "Choose the mode: 'Complete_Pipeline' for full analysis, 'Index_Calculation' for index-based analysis" + Style.RESET_ALL)

    # Parametri per la modalità 'Complete_Pipeline'
    default_group = parser.add_argument_group("Complete Pipeline Arguments")
    default_group.add_argument('-f', '--input', action="store", type=str, required=True, help=Fore.CYAN + "Specify the directory path containing the raw reads in .fastq.gz format (Required)" + Style.RESET_ALL)
    default_group.add_argument('-o', '--output', action="store", type=str, required=True, help=Fore.CYAN + "Specify the directory path where the output files will be saved (Required)" + Style.RESET_ALL)
    default_group.add_argument('-bc', '--cutpoint_based', action="store_true", default=False, help=Fore.CYAN + "Enable detection of CAG repeat numbers based on identifying the highest peaks in each of the two histograms formed by the cutpoint parameter" + Style.RESET_ALL)
    default_group.add_argument('-a', dest='amp', action='store', type=int, nargs='+', default=[5,6,7,8,9,10], help=Fore.CYAN + 'Specify amplitude values for the expected peak widths in the data (default: [5,6,7,8,9,10])' + Style.RESET_ALL)
    default_group.add_argument('-i', dest='interv', action='store', type=int, default=6, help=Fore.CYAN + 'Define the interval around identified CAG peaks for searching higher read counts (default: 6)' + Style.RESET_ALL)
    default_group.add_argument('-ti', dest='threshold_instability', action='store', type=float, default=False, help=Fore.CYAN + 'Set the relative peak height threshold for the instability index (Default False, recommended value: 0.2)' + Style.RESET_ALL)
    default_group.add_argument('-te', dest='threshold_expansion', action='store', type=float, default=False, help=Fore.CYAN + 'Set the relative peak height threshold for the Expansion index (Default False, recommended value: 0.03)' + Style.RESET_ALL)
    default_group.add_argument('-c', dest='cutpoint', action='store', type=int, default=27, help=Fore.CYAN + 'Set the cutpoint to divide the histogram into two sections (default: 27), also used to calculate AlleleRatio ' + Style.RESET_ALL)
    default_group.add_argument('-m', dest='min', action='store', type=int, default=7, help=Fore.CYAN + 'min CAG repeats (default: 7)' + Style.RESET_ALL)
    default_group.add_argument('--cag_graph', dest='cag', action='store_true', default=False, help=Fore.CYAN + 'Enable to save graphs of CAG trinucleotide repeat distributions' + Style.RESET_ALL)
    default_group.add_argument('--ccg_graph', dest='ccg', action='store_true', default=False, help=Fore.CYAN + 'Enable to save graphs of CCG trinucleotide repeat distributions' + Style.RESET_ALL)
    default_group.add_argument('--cwt', dest='cwt_finder', action='store_true', default=False, help=Fore.CYAN + 'Enable wavelet-based peak detection as an alternative to histogram-based detection' + Style.RESET_ALL)

    # Parametri per la modalità 'Index_Calculation'
    indices_group = parser.add_argument_group("Index Calculation Arguments")
    indices_group.add_argument('-p', '--path', action="store", type=str, required=False, help=Fore.CYAN + "Specify the file (in excel format) with four columns detail: Sample, CAG_Allele_1, CAG_Allele_2" + Style.RESET_ALL)

    args = parser.parse_args()

    # Logica per verificare la compatibilità dei parametri in base alla modalità
    if args.mode == "Complete_Pipeline" and not args.input:
        parser.error("In 'Complete_Pipeline' mode, you must provide the '-f' or '--input' parameter.")
    if args.mode == "Index_Calculation" and not args.path:
        parser.error("In 'Index_Calculation' mode, you must provide the '-p' or '--path' parameter with a excel file (.xlsx) containing four columns: Sample, CAG_Allele_1, CAG_Allele_2")

    # Assegnazione delle variabili
    if args.mode == "Complete_Pipeline":
        input_raw_reads = args.input + "/"
        path = args.output + "/"
        ampiezza = args.amp
        intorno = args.interv
        cag_graph = args.cag
        ccg_andWarning_graph = args.ccg
        cwt =args.cwt_finder
        cutpoint = args.cutpoint
        infMin = args.min
        biological_cutpoint = args.cutpoint_based
        ii_threshold = args.threshold_instability
        ei_threshold = args.threshold_expansion
    elif args.mode == "Index_Calculation":
        input_raw_reads = args.input + "/"
        path = args.output + "/"
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
                ii=instabilityIndex(df_distrib,cag_max_1,cag_max_2)  # instability Index , ti
                ei=expansionIndex(df_distrib,cag_max_1,cag_max_2)    # expansion index ,te
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

        df.to_excel(str(output))
    
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

            ii=instabilityIndex(df_distrib,cag_max_1,cag_max_2)  # instability Index
            ei=expansionIndex(df_distrib,cag_max_1,cag_max_2)    # expansion index

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
    
        df.to_excel(output)
    
        return df

    ################################################## RUN ###################################################################################################################################################################

    # Mostra il logo all'avvio del programma
    print_logo()

    print("Starting strmie")
    print()

    # PIPELINE COMPLETA
    if args.mode == "Complete_Pipeline":

        print("Parameters:")
        print("input directory: "+input_raw_reads)
        print("output directory: "+path)

        if biological_cutpoint:
            print("cutpoint-based: "+str(biological_cutpoint))
            print("cutpoint: "+str(cutpoint))
        elif cwt:
            print("cwt-based: "+str(cwt))
            print("width: "+str(ampiezza))
            print("interval: "+str(intorno))
        else:
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
        file_names=leggi_nomi_file_inDirectory(input_raw_reads)
        #print(file_names)

        list_data=[]
        c=0
        print("Calculate LOI and Freq.")
        for name in file_names:
            path_file=input_raw_reads+name
            if c==0:
                data=calcola_counts_and_loi(path_file,name)
                data["filename"]=name
                data=data[data.CAG_repeats>=infMin]
                c=c+1
            else:
                tmp=calcola_counts_and_loi(path_file)
                tmp["filename"]=name
                tmp=tmp[tmp.CAG_repeats>=infMin]
                data=pd.concat([data, tmp])

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
            tmp_counts.to_csv(create5+str(s)+".csv")



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
        create_html(path,create1)



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
        calculate_indices_fromFile(df_merged,campioni,path+out_indices,cutpoint=cutpoint)

        print("Done")
 
    print()
    print("The job is done, Thanks for using strmie")

if __name__ == "__main__":
    main()
