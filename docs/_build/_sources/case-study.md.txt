# Case Studies

This section provides two case studies demonstrating how to use **STRmie-HD** depending on the sequencing design:  

1. **Paired-End (PE) reads** – combine R1/R2 into single continuous reads using an external merging tool, then run STRmie-HD.  
2. **Single-End (SE) reads** – Unpaired reads (Single-End), such as Illumina or long reads from Nanopore/PacBio. Run STRmie-HD directly on gzipped FASTQ/FASTA files.  

Both flows converge to the same STRmie-HD pipeline and produce the same interactive HTML report for visual inspection.  

---

## Case Study 1 – Paired-End Reads (PE) with PEAR

When working with paired-end data, we recommend merging reads upstream with **[PEAR](https://cme.h-its.org/exelixis/web/software/pear/doc.html)** to obtain single, high-confidence sequence spanning the CAG repeat region.  

***Generic PEAR command***
```bash
pear -f forward_R1.fastq.gz -r reverse_R2.fastq.gz -v 10 -o output_prefix

```
**Parameters**
- `-f forward_R1.fastq.gz` → path to the **forward** read (R1).  
- `-r reverse_R2.fastq.gz` → path to the **reverse** read (R2).  
- `-v MIN_OVERLAP` → **minimum overlap length** (in bp) required to confidently merge a pair (default: 10).  
- `-o output_prefix` → output filename prefix; PEAR will produce files such as `output_prefix.assembled.fastq` (merged reads).


Alternative merging tool: **[FLASH (Fast Length Adjustment of Short reads)](https://github.com/ebiggers/flash)** can also be used for merging paired-end reads before running STRmie-HD.


After merging, use the merged reads as input for STRmie‑HD (see **Running the Complete Pipeline** below).

---

## Case Study 2 – Single-End Reads (SE)

For single-end datasets, no merging is required. STRmie‑HD can ingest **fastq.gz** (preferred, includes quality) or **fasta.gz** files format directly.

---


## Running the Complete Pipeline (applies to both PE-merged and SE)

Launch STRmie‑HD with default parameters on your input directory (either the PE‑merged `.assembled.fastq` outputs or the original SE files):

```bash
strmie --mode Complete_Pipeline \
       -f /path/to/input_dir \
       -o /path/to/output_dir \
       [other options]
```

**required parameters**
- `-f /path/to/input_dir/` → directory containing input reads (e.g., merged `*.assembled.fastq` from PEAR, or SE `.fastq.gz` / `.fasta.gz`).  
- `-o /path/to/output_dir/` → output directory for all STRmie‑HD outputs.


## Running the Index Calculation

> ℹ️ Guidance: If you intend to run the **Index_Calculation** mode, see the see the {ref}`html-report-step-by-step-workflow` section.
There you will find instructions on how to use the interactive HTML report to manually correct allele peaks and export a curated table for recalculating indices.


