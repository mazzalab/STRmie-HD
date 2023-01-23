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
```
Alternatively you can copy the environment from the .yml file, using the command:
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
--cag_graph Set if you want to save the CAG graphs. Default is False;
--ccg_graph Set if you want to save the CCG graphs. Default is False.
```

# Alerts
The folder containing the raw reads must contain only the input files in fastq.gz format.
