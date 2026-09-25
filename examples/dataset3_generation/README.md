# Dataset 3 generation

Scripts to reproduce the synthetic HTT exon 1 benchmark ("Dataset 3") used in
the STRmie-HD manuscript. This replaces an earlier, simpler generator that
produced fixed-length (300 bp), constant-quality reads with no genomic
flanking context; the current design embeds designed CAG/CCG/interruption
(LOI/DOI) repeat tracts into the native HTT locus sequence and simulates
realistic ONT long reads with [Badread](https://github.com/rrwick/Badread).

## 1. Build engineered references

```bash
python build_engineered_refs.py --out-dir allele_refs
```

For each of the 11 samples (`healthy_15_25_NOLOI`, ..., `very_high_expansion_18_102_LOI`),
this splices the designated CAG/CCG repeat lengths and one or more
interruption motifs (canonical, LOI, DOI, and variants) into
`htt_locus_chr4_3072000_3078000_hg38.fasta` (chr4:3072000-3078000, hg38),
replacing the endogenous repeat tract while preserving the surrounding
genomic flanking context, and writes an `allele_refs/manifest.tsv` describing
each engineered allele.

## 2. Simulate reads

```bash
# requires: pip install badread==0.4.2
./simulate_reads.sh allele_refs/manifest.tsv dataset3_reads
```

Runs Badread (v0.4.2) on each engineered reference, configured to emulate ONT
R10.4.1 chemistry (mean read length 6000 bp, sd 1500; mean identity 98%,
maximum 99.5%, sd 2; `nanopore2023` error and quality-score models; no junk,
random, or chimeric reads), then concatenates and gzips each sample's reads
into `dataset3_reads/<sample>_1.fastq.gz` (~10,000 reads/sample). These
parameters match the manuscript Methods section ("Dataset 3 synthetic
generation") exactly, so this pipeline is fully reproducible.

The resulting FASTQ files for Dataset 3 are also archived at Zenodo:
https://doi.org/10.5281/zenodo.18346811.
