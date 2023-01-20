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


parser=argparse.ArgumentParser(description="searcHD", 
                               usage="searcHD.py -f /path/to/raw_reads_file.fastq.gz -o /path/directory_output/")

parser.add_argument('-f', '--input', action="store", type=str, help="specify the directory of raw reads")
parser.add_argument('-o', '--output', action="store", type=str, help="specify the output folder")
parser.add_argument('--a', dest='amp', action='store', type=int, default=[5,6,7,8,9,10], help='This value should cover the expected width of peaks of interest')
parser.add_argument('--cag_graph', dest='cag', action='store_true', default=False, help='Set if you want to save the CAG graphs')
parser.add_argument('--ccg_graph', dest='ccg', action='store_true', default=False, help='Set if you want to save the CCG graphs')

args = parser.parse_args()
input_raw_reads=args.input
path=args.output
ampiezza=args.amp
cag_graph=args.cag
ccg_andWarning_graph=args.ccg

# In[2]:


def leggi_fastq_gz(file):
    df=pd.DataFrame(pd.read_csv(file, sep='\t', header=None,compression='gzip').values.reshape(-1, 4), columns=['ID', 'Seq', 'Sep', 'Qual'])
    return df


# In[3]:


def leggi_fastq(file):
    df=pd.DataFrame(pd.read_csv(file, sep='\t', header=None).values.reshape(-1, 4), columns=['ID', 'Seq', 'Sep', 'Qual'])
    return df


# In[4]:


def htt_exact_match(seq):
    pattern=r"(?P<CAG>(cag)+)((?P<LOI>caa)?cag){0,2}ccgcca(?P<CCG>(ccg)+)*(cct)*"

    match = re.search(pattern, seq, re.MULTILINE | re.IGNORECASE)
    if match==None:
        #print("did not find")
        return(None, None,None)
    else:
        groups = match.groups()
        groups_dict = match.groupdict()
        
        cag_len = len(match.group("CAG"))/3
        
        if match.group("CCG")==None:
            ccg_len = 0
        else:
            ccg_len = len(match.group("CCG"))/3

        is_loi = match.group("LOI") is None
        
        return(cag_len,is_loi,ccg_len)


# In[5]:


def calcola_counts_and_loi(path_file):
    
    df=leggi_fastq_gz(path_file)
    
    df_count=pd.DataFrame()
    
    cag_len=[]
    ccg_len=[]
    loi=[]
    id_sequenze=[]
    sequenze=[]

    for seq,idseq in zip(df.Seq, df.ID):    
        cag, is_loi,ccg = htt_exact_match(seq)
        cag_len.append(cag)
        ccg_len.append(ccg)
        loi.append(is_loi)
        id_sequenze.append(idseq)
        sequenze.append(seq)

    df_count["ID"]=id_sequenze
    df_count["CAG_repeats"]=cag_len
    df_count["CCG_repeats"]=ccg_len
    df_count["LOI"]=loi
    df_count["Seq"]=sequenze
    
    df_count=df_count.loc[((df_count.LOI==True) | (df_count.LOI==False))]
    
    return df_count


# In[6]:


def leggi_nomi_file_inDirectory(dir_path, estensione=False):
    import os
    res = []

    if estensione==False:
        for path in os.listdir(dir_path):
            if os.path.isfile(os.path.join(dir_path, path)):
                res.append(path)
    else :
        print("extension file: "+str(estensione))
        for path in os.listdir(dir_path):
            if path.endswith(str(estensione)):
                if os.path.isfile(os.path.join(dir_path, path)):
                    res.append(path)
    return res


# In[7]:


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
    c=0
    if len(forward)==len(reverse):
        for i,e in enumerate(forward):
            diff_letters = sum( e[j] != reverse[i][j] for j in range(len(e)) )
            if diff_letters>1:
                print("WARNING: some files do not have the same name for reverse and forward")
                c=1
    else:
        print("WARNING: some files are not paired")
        c=1
        
    if c==0: 
        print("ok, now you can see the pairs by running :")
        print()
        print("forward,reverse=coppie_PE(files)")
        print("for i in range(len(forward)):")
        print("    print(forward[i])")
        print("    print(reverse[i])")
    return forward, reverse


# In[8]:


# In[9]:


def barplot_alleli(df,titolo,name):
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize =(13, 10))
    df['CAG_repeats'].value_counts().sort_index().plot.bar()#plt.ylim(0,0.035)
    #plt.ylim(0,t[name].max()+t[name].max()/10)
    plt.rc('xtick', labelsize=7);plt.rc('ytick', labelsize=10);plt.xticks(rotation=90)
    plt.title(titolo)
    plt.xlabel("Number of CAG repetitions")
    plt.ylabel("Count of the number of repetitions")
    
    plt.savefig(name)
    plt.close(fig)


# In[10]:


def barplot_alleli_samples(df,path_outIMG):
    samples=list(df.filename.unique())
    for s in samples:
        tmp=df[df.filename==s]
        barplot_alleli(tmp,"CAG content Alleles",path_outIMG+s+".jpg")

    barplot_alleli(df,"CAG content Alleles all Samples",path_outIMG+"All_samples"+".jpg")

# In[11]:


def barplot_alleli_ccg(df,titolo,name):
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize =(13, 10))
    df['CCG_repeats'].value_counts().sort_index().plot.bar()

    plt.rc('xtick', labelsize=7);plt.rc('ytick', labelsize=10);plt.xticks(rotation=90)
    plt.title(titolo)
    plt.xlabel("Number of CCG repetitions")
    plt.ylabel("Count of the number of repetitions")
    
    plt.savefig(name)
    plt.close(fig)


# In[12]:


def instabilityIndex(df,cag1,cag2,cutOFF=True):
    
    df_pre_norm=df['CAG_repeats'].value_counts().to_frame()
    df_pre_norm=df_pre_norm.rename({'CAG_repeats':'height_peak'},axis=1)
    df_pre_norm["CAG_repeat"]=df_pre_norm.index

    if cag1==cag2:
        second_allele=df_pre_norm
    else:
        ### taglio il dataframe prendendo il punto di minimo tra i due cag, e tengo l'allele cag2
        tmp=df_pre_norm[(df_pre_norm["CAG_repeat"]>cag1)&(df_pre_norm["CAG_repeat"]<cag2)]
        soglia=tmp.height_peak.min()
        cag_soglia=tmp.CAG_repeat[tmp.height_peak==soglia].values[0]
        second_allele=df_pre_norm[df_pre_norm.CAG_repeat>cag_soglia]

    maxPeakHeight=second_allele['height_peak'][second_allele.CAG_repeat==cag2].values[0]
    sumPeakHeight=second_allele['height_peak'].sum()
    less20percent=maxPeakHeight*0.20
    if cutOFF:
        df_pre_norm=second_allele[second_allele['height_peak']>less20percent]
    df_norm=df_pre_norm[['height_peak']].div(df_pre_norm[['height_peak']].sum(axis=0), axis=1)
    df_norm["CAG_repeat"]=df_pre_norm["CAG_repeat"]
    
    #max_repeatition=df_norm["CAG_repeat"].iloc[0]
    pre=len(df_norm[df_norm["CAG_repeat"]<cag2])#max_repeatition --> cag2
    post=len(df_norm[df_norm["CAG_repeat"]>cag2])#max_repeatition --> cag2
    pre_ord= [*range(-pre, 0, 1)]
    post_ord= [*range(0, post+1, 1)]
    ordinamento=pre_ord+post_ord
    df_norm=df_norm.sort_values('CAG_repeat')
    df_norm["ordinamento"]=ordinamento
    df_norm["heightXchanges"]=df_norm['height_peak']*df_norm["ordinamento"]

    instability_index=df_norm["heightXchanges"].sum()
    #print(f'instability_index: {instability_index}')

    return df_norm


# In[13]:


def expansionIndex(df,cag1,cag2,cutOFF=True):
    
    df_pre_norm=df['CAG_repeats'].value_counts().to_frame()
    df_pre_norm=df_pre_norm.rename({'CAG_repeats':'height_peak'},axis=1)
    df_pre_norm["CAG_repeat"]=df_pre_norm.index

    if cag1==cag2:
        second_allele=df_pre_norm
    else:
        ### taglio il dataframe prendendo il punto di minimo tra i due cag, e tengo l'allele cag2
        tmp=df_pre_norm[(df_pre_norm["CAG_repeat"]>cag1)&(df_pre_norm["CAG_repeat"]<cag2)]
        soglia=tmp.height_peak.min()
        cag_soglia=tmp.CAG_repeat[tmp.height_peak==soglia].values[0]
        second_allele=df_pre_norm[df_pre_norm.CAG_repeat>cag_soglia]
  
    maxPeakHeight=second_allele['height_peak'][second_allele.CAG_repeat==cag2].values[0]
    sumPeakHeight=second_allele['height_peak'].sum()
    
    less3percent=maxPeakHeight*0.03
    
    indice_maxPeakHeight=second_allele[second_allele["height_peak"]==maxPeakHeight]["CAG_repeat"].values[0]
    df_pre_norm=second_allele[second_allele['CAG_repeat']>=indice_maxPeakHeight]
    
    if cutOFF:
        df_pre_norm=df_pre_norm[df_pre_norm["height_peak"]>less3percent]
    
    df_norm=df_pre_norm[['height_peak']].div(maxPeakHeight)
    df_norm["CAG_repeat"]=df_pre_norm["CAG_repeat"]
        
    ordinamento=[*range(0, len(df_norm), 1)]
    
    df_norm=df_norm.sort_values('CAG_repeat')
    df_norm["ordinamento"]=ordinamento
    df_norm["heightXchanges"]=df_norm['height_peak']*df_norm["ordinamento"]
    
    expansion_index=df_norm["heightXchanges"].sum()
    
    return df_norm


# In[14]:


def make_report_indices(df_II,df_EI,sample):
    ei=df_EI.heightXchanges.sum()
    ii=df_II.heightXchanges.sum()
    
    # Value of the max number of CAG repetitions: 15.0

    cag_maxRepeats=df_EI["CAG_repeat"][df_EI.heightXchanges==0].values[0]
    check=df_II["CAG_repeat"][df_II.heightXchanges==0].values[0] # controllo che siano uguali
    
    if cag_maxRepeats!=check: print("Ambiguity in the maximum peak of sample:"+sample+":"+str(cag_maxRepeats)+","+str(check))
    
    ii_height_peak_norm=df_II["height_peak"][df_II.heightXchanges==0].values[0]
    
    return ii, ei, cag_maxRepeats, ii_height_peak_norm


# In[15]:


def find_peaks_two_alleles(df, ampiezza=6, intorno=6):
    from scipy import signal

    tmp=df['CAG_repeats'].value_counts().sort_index().to_frame()
    
    #### (1) filtro tutte le ripetizioni CAG minori di 9
    tmp=tmp[tmp.index>8]

    if tmp.empty:
        print("ERROR_1")
        print("The coverage of the sample: "+str(df.filename.values[0])+" is not sufficient to perform the analysis.") 
        print("Remove it from the folder and try to run searcHD on the remaining samples.")
        print("#######")

    altezze=tmp.CAG_repeats.values
    max_peak=int(list(tmp[tmp.CAG_repeats==tmp.CAG_repeats.max()].index)[0])
    peak_indices= signal.find_peaks_cwt(altezze,widths=[5,6,7,8,9,10])

    if len(peak_indices)==2: # controllo se i picchi trovati sono 2
        peaks=[int(tmp[tmp.CAG_repeats==altezze[peak_indices[0]]].index[0]),int(tmp[tmp.CAG_repeats==altezze[peak_indices[1]]].index[0])]
        ## controllo il primo picco
        if (peaks[0]==max_peak) & (peaks[0]<36) & (peaks[1]>35): # 26
            
            # controllo il secondo picco se è il massimo in un intorno e prendo il massimo
            df_check_peak_2=df.loc[(df.CAG_repeats<(peaks[1]+intorno)) & (df.CAG_repeats>(peaks[1]-intorno))]
            check_cag_2=list(df_check_peak_2.CAG_repeats.value_counts().index)[0]
            
            return peaks[0],check_cag_2
        
        else:
            df_check_peak_1=df.loc[(df.CAG_repeats<(peaks[0]+intorno)) & (df.CAG_repeats>(peaks[0]-intorno))]
            check_cag_1=list(df_check_peak_1.CAG_repeats.value_counts().index)[0]
            ## controllo il primo picco
            if (max_peak == check_cag_1) & (peaks[0]<36) & (peaks[1]>35): #26
                
                # controllo il secondo picco se è il massimo
                df_check_peak_2=df.loc[(df.CAG_repeats<(peaks[1]+intorno)) & (df.CAG_repeats>(peaks[1]-intorno))]
                check_cag_2=list(df_check_peak_2.CAG_repeats.value_counts().index)[0]

                return check_cag_1,check_cag_2
            else:
                return "warning","warning"
  
    else:
        return "warning","warning"


# In[16]:


def report_to_excel(data,campioni,output,ampiezza=6):
    instInd=[]
    expInd=[]
    cag=[]
    iiHeight=[]
    cag_max_alleles_1=[]
    cag_max_alleles_2=[]
    percentage=[]
    
    for c in campioni:
        data_campione=data[data.filename==c]
        
        cag_max_1, cag_max_2 = find_peaks_two_alleles(data_campione, ampiezza)  # picchi massimi nei due alleli binomial distribution
        cag_max_alleles_1.append(cag_max_1)
        cag_max_alleles_2.append(cag_max_2)
        
        if type(cag_max_1)==str:
            instInd.append("warning")
            expInd.append("warning")
            cag.append("warning")
            iiHeight.append("warning")
            percentage.append("warning")

        else:
            t1=instabilityIndex(data_campione,cag_max_1,cag_max_2)  # instability Index
            t2=expansionIndex(data_campione,cag_max_1,cag_max_2)    # expansion index
            ii, ei, cag_maxRepeats, ii_height_peak_norm=make_report_indices(t1,t2,c)
            instInd.append(ii)
            expInd.append(ei)
            cag.append(cag_maxRepeats)
            iiHeight.append(ii_height_peak_norm)
        
            ## LOI: percentage loi sulle reads espanse dell'allele mutato
            df_loi=data_campione[data_campione.CAG_repeats>=cag_max_2]
            if ((True in df_loi.LOI.values) & (False in df_loi.LOI.values)):         # solo alcuni hanno la loi
                loi=df_loi.LOI.value_counts().to_frame().loc[True][0]
                nonLoi=df_loi.LOI.value_counts().to_frame().loc[False][0]
                percentage.append((loi/(loi+nonLoi))*100)
            elif (True in df_loi.LOI.values) & ((False in df_loi.LOI.values)==False): # hanno tutti la loi
                percentage.append(100)
            else:                                                                      # nessuno ha la loi
                percentage.append(0)

        
    df=pd.DataFrame()
    df["Sample"]=campioni
    df["Instability_Index"]=instInd
    df["Expansion_Index"]=expInd
    df["CAG_repeatsPeak"]=cag
    df["LOI_percentage"]=percentage
    df["InstabilityPeakNormalized"]=iiHeight
    df["CAG_repeatsPeak_Allele_1"]=cag_max_alleles_1
    df["CAG_repeatsPeak_Allele_2"]=cag_max_alleles_2
    
    df.to_excel(str(output))
    
    return df


# In[17]:


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


# In[18]:


def force_search(t): # t corrisponde al data_campione presente nella funzione report_to_excel
    name=t.filename

    t=t.CAG_repeats.value_counts().to_frame()
    t["CAG_repeat"]=t.index
    t=t.rename({'CAG_repeats':'height_peak'},axis=1)
    t=t[t.CAG_repeat>3]
    t.sort_values(by=["CAG_repeat"],inplace=True)
    media=t.height_peak.mean()
    massimo=t.height_peak.max()
    cag_max=t.CAG_repeat[t.height_peak==massimo].values[0]
    
    c=3
    if cag_max<=27: # mi muovo verso destra
        for i in list(t.CAG_repeat[t.CAG_repeat>cag_max])[3:]: ### skippo i primi tre picchi vicini
            
            picco_dx=t.height_peak[t.CAG_repeat==i].values[0]
            if picco_dx<int(media):
                new_cag=t.CAG_repeat[t.height_peak==picco_dx].values[0]
                t_new=t[t.CAG_repeat>new_cag]
                second_massimo=t_new.height_peak.max()
                ### eccezione
                if math.isnan(second_massimo):
                    altezze=t.CAG_repeat.values
                    peak_indices= signal.find_peaks_cwt(altezze,widths=6)
                    second_cag_max=int(t[t.CAG_repeat==altezze[peak_indices[0]]].index[0])
                    break
                ###
                second_cag_max=t_new.CAG_repeat[t_new.height_peak==second_massimo].values[0]
                break

    elif cag_max<=36: # controllo a destra e sinistra
        for i in list(t.CAG_repeat[t.CAG_repeat>cag_max])[3:]: ### destra , ### skippo i primi tre picchi vicini
            picco_dx=t.height_peak[t.CAG_repeat==i].values[0]
            if picco_dx<int(media):
                p1_new_cag=t.CAG_repeat[t.height_peak==picco_dx].values[0]
                p1_t_new=t[t.CAG_repeat>p1_new_cag]
                p1_second_massimo=p1_t_new.height_peak.max()
                ### eccezione
                if math.isnan(p1_second_massimo):
                    altezze=t.CAG_repeat.values
                    peak_indices= signal.find_peaks_cwt(altezze,widths=6)
                    p1_second_cag_max=int(t[t.CAG_repeat==altezze[peak_indices[0]]].index[0])
                    break
                ###
                p1_second_cag_max=p1_t_new.CAG_repeat[p1_t_new.height_peak==p1_second_massimo].values[0]

        c=3
        for i in list(t.CAG_repeat[t.CAG_repeat<cag_max])[:-3]: ### sinistra, ### skippo i primi tre picchi vicini
            picco_sx=t.height_peak[t.CAG_repeat==i].values[0]
            if picco_sx<int(media):
                p2_new_cag=t.CAG_repeat[t.height_peak==picco_sx].values[0]
                p2_t_new=t[t.CAG_repeat<p2_new_cag]
                p2_second_massimo=p2_t_new.height_peak.max()
                ### eccezione
                if math.isnan(p2_second_massimo):
                    altezze=t.CAG_repeat.values
                    peak_indices= signal.find_peaks_cwt(altezze,widths=6)
                    p2_second_cag_max=int(t[t.CAG_repeat==altezze[peak_indices[0]]].index[0])
                    break
                ###
                p2_second_cag_max=p2_t_new.CAG_repeat[p2_t_new.height_peak==p2_second_massimo].values[0]

        if p1_second_massimo> p2_second_massimo:
            second_massimo=p1_second_massimo
            second_cag_max=p1_second_cag_max
        else: 
            second_massimo=p2_second_massimo
            second_cag_max=p2_second_cag_max

    else:  # mi muovo verso sinistra
        for i in list(t.CAG_repeat[t.CAG_repeat<cag_max])[:-3]: ### skippo i primi tre picchi vicini
            picco_sx=t.height_peak[t.CAG_repeat==i].values[0]
            if picco_sx<int(media):
                new_cag=t.CAG_repeat[t.height_peak==picco_sx].values[0]
                t_new=t[t.CAG_repeat<new_cag] ### non è detto che esista, potrebbe essere l'estremo del df
                second_massimo=t_new.height_peak.max()
                ### eccezione
                if math.isnan(second_massimo):
                    print(name)
                    altezze=t.CAG_repeat.values
                    peak_indices= signal.find_peaks_cwt(altezze,widths=6)
                    second_cag_max=int(t[t.CAG_repeat==altezze[peak_indices[0]]].index[0])
                    print(second_cag_max)
                    break
                ###
                second_cag_max=t_new.CAG_repeat[t_new.height_peak==second_massimo].values[0]
                break
    
    try: 
        tmp=second_cag_max
    except NameError:
        second_cag_max= cag_max
        
    if cag_max>second_cag_max:
        cag1=second_cag_max
        cag2=cag_max
    else:
        cag1=cag_max
        cag2=second_cag_max
    return cag1,cag2


# In[19]:


def force_findingPeaks(data,campioni,output):
    instInd=[]
    expInd=[]
    cag=[]
    iiHeight=[]
    cag_max_alleles_1=[]
    cag_max_alleles_2=[]
    percentage=[]
    
    for c in campioni:
        data_campione=data[data.filename==c]
        
        cag_max_1, cag_max_2 = force_search(data_campione)
        cag_max_alleles_1.append(cag_max_1)
        cag_max_alleles_2.append(cag_max_2)
        
        t1=instabilityIndex(data_campione,cag_max_1,cag_max_2)  # instability Index
        t2=expansionIndex(data_campione,cag_max_1,cag_max_2)    # expansion index

        ii, ei, cag_maxRepeats, ii_height_peak_norm=make_report_indices(t1,t2,c)
        instInd.append(ii)
        expInd.append(ei)
        cag.append(cag_maxRepeats)
        iiHeight.append(ii_height_peak_norm)
        
        ## LOI: percentage loi sulle reads espanse dell'allele mutato
        df_loi=data_campione[data_campione.CAG_repeats>=cag_max_2]
        if ((True in df_loi.LOI.values) & (False in df_loi.LOI.values)):         # solo alcuni hanno la loi
            loi=df_loi.LOI.value_counts().to_frame().loc[True][0]
            nonLoi=df_loi.LOI.value_counts().to_frame().loc[False][0]
            percentage.append((loi/(loi+nonLoi))*100)
        elif (True in df_loi.LOI.values) & ((False in df_loi.LOI.values)==False): # hanno tutti la loi
            percentage.append(100)
        else:                                                                      # nessuno ha la loi
            percentage.append(0)
        
    df=pd.DataFrame()
    df["Sample"]=campioni
    df["Instability_Index"]=instInd
    df["Expansion_Index"]=expInd
    df["CAG_repeatsPeak"]=cag
    df["LOI_percentage"]=percentage
    df["InstabilityPeakNormalized"]=iiHeight
    df["CAG_repeatsPeak_Allele_1"]=cag_max_alleles_1
    df["CAG_repeatsPeak_Allele_2"]=cag_max_alleles_2
    
    df.to_excel(output)
    
    return df




################################################## RUN
print("Starting searcHD")


# In[25]:


dir1="CAG_graphs"
dir2="CCG_alleles_graphs"
dir3="warning_case"
dir4="forced_graphs"

create1 = os.path.join(path, dir1)
create2 = os.path.join(path, dir2)
create3 = os.path.join(path, dir3)
create4 = os.path.join(path, dir4)
createFolder=True

if cag_graph & ccg_andWarning_graph:
    folders=[create1,create2,create3,create4]
elif cag_graph:
    folders=[create1]
elif ccg_andWarning_graph:
    folders=[create2,create3,create4]
else:
    createFolder=False

if createFolder:
    for c in folders:
        try:
            os.mkdir(c)
            print("Directory '%s' created" % c)
        except FileExistsError:
            print("Directory '%s' already exists" % c)

# In[26]:


create1=create1+"/"
create2=create2+"/"
create3=create3+"/"
create4=create4+"/"


# In[27]:

print("create dataframe from raw reads")
file_names=leggi_nomi_file_inDirectory(input_raw_reads,".fastq.gz")

c=0
for name in file_names:
    path_file=input_raw_reads+name
    if c==0:
        data=calcola_counts_and_loi(path_file)
        data["filename"]=name
        c=c+1
    else:
        tmp=calcola_counts_and_loi(path_file)
        tmp["filename"]=name
        data=pd.concat([data, tmp])

data["len"]=data.Seq.apply(lambda x: len(x))

if cag_graph:
    print("plotting the cag graphs")
    barplot_alleli_samples(data,create1)

campioni=list(data.filename.unique())

out1="report_0.xlsx"
print("calculate cag content, indices and draft report")
pear_dataframe=report_to_excel(data,campioni,path+out1,ampiezza)

print("calculate ccg content")
final,reRun=ccg_count(data,pear_dataframe,create2,create3,ccg_andWarning_graph)


if reRun.empty:
    outFile="Final_report.xlsx"
    print("save final report")
    final.to_excel(path+outFile,index=False)
    print("The job is done, Thanks for using searcHD")
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
    print("save final report")
    final.to_excel(path+outFile,index=False)
    print("The job is done, Thanks for using searcHD")

