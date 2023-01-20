# searcHD
Automated Huntington Disease polyQ pattern scanner

# Installing searcHD
```
- Prerequisites: Python version >= Python 3.10.4

- Install dependency using conda
    conda install -c conda-forge regex
    conda install -c anaconda numpy
    conda install -c anaconda pandas 
    conda install -c conda-forge matplotlib
    conda install mathematical
    conda install -c anaconda scipy (alternative channel: conda-forge)
```
# Usage
```
python searcHD.py -f /path/directory/files/.fastq.gz -o /path/directory/output/
```
# Optional 

```
--a (integer or list of integer). This value should cover the expected width of peaks of interest. Default is a list [5,6,7,8,9,10];
--cag_graph Set if you want to save the CAG graphs. Default is False.
--ccg_graph Set if you want to save the CCG graphs. Default is False.
```

# Alerts
The folder containing the raw reads must contain only the input files in fastq.gz format.
