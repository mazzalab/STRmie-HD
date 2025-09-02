# Input file types

STRmie-HD is optimized for **high-quality targeted sequencing data** (e.g., Illumina MiSeq, PacBio SMRT).  
Its performance may be reduced on **lower-coverage or noisier platforms** (such as Oxford Nanopore) or in **whole-genome/low-depth clinical pipelines**, where additional parameter tuning or error-correction strategies may be required.  

Because the tool relies on **exact regular-expression matching** to identify repeat tracts and interruption motifs, reads with sequencing errors or incomplete coverage of the repeat region may be excluded from the analysis.  
As a result, STRmie-HD performs best when applied to datasets with **high sequencing accuracy and sufficient coverage** of the target locus.

STRmie-HD accepts as input:  
- **fastq.gz** → gzipped file of raw reads with base calls and quality scores.  
- **fasta.gz** → gzipped file of sequences without quality scores.  

> 💡 *Tip*: Merging is especially useful for repeat-rich regions (e.g., *HTT*), where overlap-aware merging improves accuracy by using base qualities and positional evidence.  
