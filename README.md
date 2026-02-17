# STRmie-HD

[![Bioconda](https://anaconda.org/bioconda/strmie-hd/badges/version.svg)](https://anaconda.org/bioconda/strmie-hd)
[![Bioconda downloads](https://anaconda.org/bioconda/strmie-hd/badges/downloads.svg)](https://anaconda.org/bioconda/strmie-hd)
[![BioContainer (Quay)](https://img.shields.io/badge/quay.io-biocontainers%2Fstrmie--hd-blue)](https://quay.io/repository/biocontainers/strmie-hd)
[![Docs](https://img.shields.io/badge/docs-GitHub%20Pages-brightgreen)](https://mazzalab.github.io/STRmie-HD/)


**STRmie-HD** (Short Tandem Repeat Mapping and Identification Engine – Huntington's Disease) is an interactive, Python-based tool designed to support the curation, visualization, and interpretation of short tandem repeat (STR) genotyping data obtained from Huntington's Disease (HD) patients. It enables the prediction, refinement and validation of CAG/CCG repeat expansion results in the context of HD, by highlighting cases of allelic instability or potential misclassification.

## Key Features
- Analysis and curation of STR Huntington's Disease genotypes (CAG/CCG)
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


### 🔹 2. Index Calculation Only

This mode is intended for situations where **automatic allele calling requires manual adjustment**.  
When inspecting the interactive HTML report, the user may decide that the automatically identified allele peaks are not accurate.  
Through the HTML interface, it is possible to:
- Visually explore **CAG and CCG histograms** for each sample.  
- **Manually adjust allele peak values** if necessary.  
- Export a curated Excel table with the corrected allele definitions.  

The ***Index_Calculation*** mode takes this adjusted table as input and **recomputes all instability and expansion indices (II, EI)** accordingly.  
This ensures that downstream results are based on manually validated allele assignments.

```bash
strmie --mode Index_Calculation \
       -f /path/to/input_dir \
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

> *[Paper title here]*  
> [Alessandro Napoli, Niccolò Liorni, Tommaso Mazza]  
> *Journal*, [Year], DOI: [insert DOI]

