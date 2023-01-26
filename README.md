# searcHD
Automated Huntington Disease polyQ pattern scanner

# Installing searcHD
```
- Prerequisites: Python version >= 3.10

- Install dependency:
    conda create --name searcHD_env python>=3.10
    conda activate searcHD_env
    conda install -c conda-forge regex
    conda install -c anaconda numpy
    conda install -c anaconda pandas 
    conda install -c conda-forge matplotlib
    pip install python-math
    conda install -c conda-forge scipy
    conda install -c conda-forge libstdcxx-ng
    conda install -c anaconda openpyxl 
```
Alternatively you can download the searcHD_env.yml file on github, and copy the environment using the command:
```
conda env create -n searcHD_env -f searcHD_env.yaml
```
# Usage
```
python searcHD.py -f /path/directory/files/.fastq.gz -o /path/directory/output/
```
# Optional 

```
--a (integer or list of integer). This value should cover the expected width of peaks of interest. Default is a list [5,6,7,8,9,10];
--cag_graph Set if you want to save the CAG graphs;
--ccg_graph Set if you want to save the CCG graphs.
```

# Alerts
The folder containing the raw reads must contain only the input files in fastq.gz format.

# Paired-End sequencing
The use of Paired-End sequencing provides an opportunity to generate much longer reads by overlapping and merging read pairs before using searcHD. Furthermore this would guarantee only one output for each sample instead of two (forward and reverse). 
Of course, you need to make sure that your sequencing is designed for minimal ovarlapping between forward and reverse reads.
In this case it may be useful to use tools for merging PE reads. Let's see an example using the PEAR tool:
```
pear  -f {forward_R1.fastq.gz} -r {reverse_R2.fastq.gz} -v {min_overlap} -o {path_output}
```
PEAR produces an output file with extension ".assembled.fastq" which can be used (after compressing it into fastq.gz) as an input file for searcHD.
