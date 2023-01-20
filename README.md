# searchHD
Automated Huntington Disease polyQ pattern scanner

# usage
python searcHD.py -f /path/directory/files/.fastq.gz -o /path/directory/output/

# optional 
--a (integer or list of integer). This value should cover the expected width of peaks of interest. Default is a list [5,6,7,8,9,10]
--cag_graph Set if you want to save the CAG graphs. Default is False.
--ccg_graph Set if you want to save the CCG graphs. Default is False.

# Alerts
The folder containing the raw reads must contain only the input files in fastq.gz format.
