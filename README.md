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
- `strmie-repeat`: a separate, catalog-driven submodule extending the same alignment-free counting approach to other repeat-expansion loci (FMR1, C9orf72; user-definable via JSON), see below

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

### 🔹 3. Extending Beyond HTT: `strmie-repeat`

STRmie-HD's core idea — counting a tandem repeat motif directly in raw reads, with no alignment step
required — is not inherently HTT-specific. `strmie-repeat` is a separate, additive CLI entry point that
applies this same alignment-free approach to other repeat-expansion loci, without touching or risking any
regression in the validated HTT pipeline (`strmie`/`strmie --mode ...`): it shares no code path with
`strmie/main.py`, reusing only two already-generic utility functions (peak calling, reverse-complement).

Two loci ship built in:

| Locus | Disease | Repeat motif | Pattern |
|---|---|---|---|
| `FMR1` | Fragile X Syndrome | CGG | Curated, interruption-aware: sums CGG repeat units across up to two AGG interruptions, the same way `strmie`'s own HTT engine is built from HTT's real repeat architecture rather than a generic motif-repeat. Loss of these AGG interruptions is itself the established marker of CGG instability/expansion risk, directly analogous to HTT's own loss-of-interruption (LOI) concept. |
| `C9orf72` | ALS/FTD | GGGGCC | Plain motif-repeat match. Unlike FMR1/HTT, C9orf72 has no comparably well-characterized internal interruption architecture in the literature to curate a compound pattern from — pathogenic alleles are essentially one long uninterrupted GGGGCC run, matching how repeat-primed PCR sizes this locus clinically. |

Both patterns were derived and verified against real GRCh38 reference sequence (see
`strmie/scripts/repeat_catalog.py` for the coordinates and derivation notes), not written generically.

```bash
strmie-repeat -f /path/to/sample.fastq.gz --locus FMR1 -o /path/to/output_dir
strmie-repeat -f /path/to/sample.fastq.gz --locus C9orf72 -o /path/to/output_dir
```

`-f`/`--input` accepts one or more `.fastq.gz`/`.fasta.gz` files and/or directories, same as `strmie`.
Output mirrors the canonical HTT pipeline's, minus what's HTT-specific (no CCG tract, no LOI/DOI
interruption tracking): an Excel report (`Generic_repeat_report.xlsx`, one row per sample with the two
called alleles plus the somatic **Instability Index (II)** and **Expansion Index (EI)**, computed the same
way as `strmie`'s own indices), a self-contained interactive `report.html` (cohort table, per-sample
histogram, glossary — same visual style as the HTT pipeline's report), and a per-sample
`<sample>.repeat_counts.csv` of every read's repeat count. `-ti`/`-te` set the same PCR-noise filtering
thresholds as `strmie`'s own `-ti`/`-te` flags.

**Adding your own locus:** `--locus` also accepts a path to a JSON file instead of a built-in name, with at
minimum a `gene` and `motif` field (a plain motif-repeat match, same as C9orf72 above). Optional fields:
`min_repeat_units`, `reference_build`/`region` (reused as-is by `--bam`/`--cram` extraction), and, for a
locus with its own known interruption architecture, `pattern`/`repeat_unit_groups` for a curated compound
regex (see the FMR1 entry in `repeat_catalog.py` for the schema).

**Noisy long reads:** add `--nanopore` for the same fuzzy, error-tolerant matching design as `strmie`'s own
`--nanopore` flag (see the Complete Pipeline section above), fuzzy-matching short flanking sequences with a
bounded edit distance and estimating repeat length from the span between them, rather than requiring an
unbroken run of the motif. This is a separate, fresh implementation for arbitrary catalog loci (not a reuse
of `strmie`'s HTT-specific nanopore code) and only works for a locus that defines `flank_upstream`/
`flank_downstream` (built-in FMR1/C9orf72 do). Default (no `--nanopore`) remains exact matching, unaffected
by this flag. `--np-max-roi`/`--np-max-edits`/`--np-seed-len`/`--np-no-bestmatch` tune it (mirroring
`strmie`'s own `--np-*` names); raise `--np-max-roi` for very large expansions, and consider
`--np-no-bestmatch` for large cohorts, since `regex`'s optimal-alignment search can get slow at wide ROI/edit
settings on many reads.

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
