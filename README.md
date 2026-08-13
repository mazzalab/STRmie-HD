<p align="center">
  <img
    src="docs/_images/bfx_logo.png"
    alt="Fondazione LIRH logo"
    height="150"
  >
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
  <img
    src="docs/_images/logo lirh RESTYLING pos.png"
    alt="BFX logo"
    height="150"
  >
</p>

<p align="center">
  Developed in collaboration with
  <strong>Fondazione LIRH – Lega Italiana Ricerca Huntington</strong>.
</p>

# STRmie-HD

[![Bioconda](https://anaconda.org/bioconda/strmie-hd/badges/version.svg)](https://anaconda.org/bioconda/strmie-hd)
[![Bioconda downloads](https://anaconda.org/bioconda/strmie-hd/badges/downloads.svg)](https://anaconda.org/bioconda/strmie-hd)
[![BioContainer (Quay)](https://img.shields.io/badge/quay.io-biocontainers%2Fstrmie--hd-blue)](https://quay.io/repository/biocontainers/strmie-hd)
[![Docs](https://img.shields.io/badge/docs-GitHub%20Pages-brightgreen)](https://mazzalab.github.io/STRmie-HD/)


**STRmie-HD** (Short Tandem Repeat Mapping and Identification Engine – Huntington's Disease) is an interactive, Python-based tool designed to support the curation, visualization, and interpretation of short tandem repeat (STR) genotyping data obtained from Huntington's Disease (HD) patients. It enables the prediction, refinement and validation of CAG/CCG repeat expansion results in the context of HD, by highlighting cases of allelic instability or potential misclassification.

**Supported data:** short-read (Illumina) and long-read (PacBio HiFi/CCS, Oxford Nanopore) sequencing platforms, taking either raw reads (`.fastq.gz`/`.fasta.gz`) or, if you only have aligned data, indexed BAM/CRAM directly. Nanopore's higher per-read error rate needs the `--nanopore` flag (see the Complete Pipeline section below); by default STRmie-HD uses exact matching, best suited to Illumina/PacBio-level accuracy.

## Key Features
- Analysis and curation of STR Huntington's Disease genotypes (CAG/CCG)
- Multi-platform: Illumina, PacBio, Oxford Nanopore; FASTQ/FASTA or BAM/CRAM input
- Calculation of Instability Index (II) and Expansion Index (EI)
- Local graphical interface for manual inspection (HTML report)

## Documentation
The full documentation is available here:  
👉 [STRmie-HD Documentation](https://mazzalab.github.io/STRmie-HD/)

---

## 🛠️ Installation

### 1. Clone the repository

```bash
git clone https://github.com/mazzalab/STRmie-HD.git
cd STRmie-HD
```

### 2. Create and activate the conda environment

```bash
conda env create -f STRmie.yml
conda activate STRmie
```

### 3. Install the package in development mode

```bash
pip install -e .
```

### 4. Run the CLI

```bash
strmie --help
```

---

## ✅ Automated Testing with Pytest

The project includes a test suite that validates the core functionalities of both operational modes using example input and expected output files.  
This ensures the tool works as intended after installation or modification.

### 🔸 Run tests with:

```bash
pytest tests/test_strmie.py
```

---


## 🧪 Command-line Usage

STRmie-HD provides two main operational modes: **Complete_Pipeline** and **Index_Calculation**. 

---

### 🔹 1. Complete Pipeline

This mode executes the **entire workflow** starting directly from raw sequencing data (`.fastq.gz` or `.fasta.gz`).  
It automatically performs:
- Histogram-based **CAG allele peak calling**.  
- **CCG repeat assignment** for each allele.  
- Calculation of **Instability and Expansion indices (II, EI)**.  
- Detection of **Loss of Interruption (LOI)** and **Duplication of Interruption (DOI)** events.  
- Generation of both an **Excel summary report** and an **interactive HTML report** for visual inspection.

```bash
strmie --mode Complete_Pipeline \
       -f /path/to/input_dir \
       -o /path/to/output_dir \
       [other options]
```

`-f`/`--input` also accepts a single file, or a space-separated list of files and/or directories, instead of a whole directory, to process only a subset of samples, e.g. `-f /path/to/sample1.fastq.gz /path/to/sample2.fastq.gz`.

**Paired-end data:** R1/R2 must be merged into a single sequence per sample before analysis; passing
both files directly processes each as its own (incorrect) sample. Add `--merge_paired_end` to detect
R1/R2 pairs by filename (`*_R1`/`*_R2` or `*_1`/`*_2`) and merge each pair with [PEAR](https://anaconda.org/bioconda/pear)
automatically (requires `pear` on PATH: `conda install -c bioconda pear`). Files not part of a detected
pair are still processed individually, and the flag is off by default.

**BAM/CRAM input:** if your reads are already aligned and you don't have (or don't want to regenerate)
FASTQ files, use `--bam`/`--cram` instead of `-f`/`--input` — these are separate, dedicated flags, not
extensions accepted by `-f`, and exactly one of `-f`/`--bam`/`--cram` must be given. Each accepts one or
more indexed BAM/CRAM files or directories, same as `-f`. Reads near the HTT locus are extracted
internally (including reads that failed to align well, which is common for large repeat expansions) and
converted to FASTQ before running the same pipeline unchanged — no need for `--merge_paired_end`, paired
reads extracted this way are merged automatically.

```bash
strmie --mode Complete_Pipeline \
       --bam /path/to/sample.bam \
       -o /path/to/output_dir

# CRAM requires the reference it was aligned against, to decode it
strmie --mode Complete_Pipeline \
       --cram /path/to/sample.cram --reference /path/to/reference.fa \
       -o /path/to/output_dir
```

The reference build (GRCh38 or GRCh37/hg19) is auto-detected from the file header; override with
`--locus-build` or, for a non-HTT locus or non-standard reference, `--bam-region chrom:start-end`.
Requires `pysam` (`pip install strmie[bam]`, or `conda install -c bioconda pysam`) and `samtools` on
PATH (`conda install -c bioconda samtools`).


### 🔹 2. Index Calculation Only

This mode is intended for situations where **automatic allele calling requires manual adjustment**.  
When inspecting the interactive HTML report, the user may decide that the automatically identified allele peaks are not accurate.  
Through the HTML interface, it is possible to:
- Visually explore **CAG and CCG histograms** for each sample.  
- **Manually adjust allele peak values** if necessary.  
- Export a curated Excel table with the corrected allele definitions.  

The ***Index_Calculation*** mode takes this adjusted table as input and **recomputes all instability and expansion indices (II, EI)** accordingly.  
This ensures that downstream results are based on manually validated allele assignments.

`-o` must point at **the same output directory** used for the original Complete Pipeline run, since this mode
reuses the raw per-read counts already saved in its `raw_counts` subfolder. `-f`/`--input` is not needed in this mode.

```bash
strmie --mode Index_Calculation \
       -o /path/to/output_dir \
       -p /path/to/CAG_data_for_recalculating_indices.xlsx
```

---


## 📦 Dependencies

Included in the `STRmie.yml` file. Core packages include:

- `pandas`, `numpy`, `openpyxl`, `xlsxwriter`
- `pytest` (for running the test suite)
- JavaScript frontend tools: `Chart.js`, `Bootstrap`, `DataTables`, `XLSX.js`

---

## 📄 License

This project is licensed under the MIT License. See `LICENSE` for details.

---

## 👩‍🔬 Authors

Developed by the [Mazzalab](https://github.com/mazzalab), Italy.  
For questions, please open an issue or contact the maintainers.

---

## 🔗 Citation

If you use STRmie-HD, please cite:

> *[STRmie-HD: A Short Tandem Repeat Mapping and Identification Engine for Interruption-Aware Genotyping and Somatic Mosaicism Profiling in Huntington’s Disease]*  
