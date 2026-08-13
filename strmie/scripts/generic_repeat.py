#!/usr/bin/env python
# coding: utf-8
"""
Generic, catalog-driven repeat counting -- STRmie-HD's core alignment-free
matching idea (find a repeated k-mer motif directly in raw reads, no
alignment needed), generalized to an arbitrary locus/motif rather than the
hardcoded HTT CAG/CCG regex in strmie/scripts/pattern.py.

This module is intentionally self-contained: it does not import from or
modify strmie/main.py or the HTT-specific matching/report code, so it
carries zero risk to the existing, validated HD pipeline. It reuses only
strmie.scripts.peaks._two_highest_peaks (the topographic-prominence peak
picker), which is already parametrized over an arbitrary count Series and
unmodified here, and strmie.scripts.pattern.get_fast_rev_comp, a small,
also-unmodified utility.
"""

import gzip
import os
import re

import pandas as pd

from strmie.scripts.peaks import _two_highest_peaks
from strmie.scripts.pattern import get_fast_rev_comp
from strmie.scripts.repeat_catalog import load_catalog
from strmie.scripts.generic_indices import somatic_indices


def _read_fastq_or_fasta(path):
    """Yield (read_id, sequence) for a .fastq.gz or .fasta.gz file."""
    is_fastq = ".fastq" in path
    with gzip.open(path, "rt") as f:
        if is_fastq:
            while True:
                header = f.readline()
                if not header:
                    break
                seq = f.readline().strip()
                f.readline()  # '+'
                f.readline()  # quality
                yield header.strip().lstrip("@"), seq
        else:
            read_id, seq_chunks = None, []
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if read_id is not None:
                        yield read_id, "".join(seq_chunks)
                    read_id, seq_chunks = line.lstrip(">"), []
                else:
                    seq_chunks.append(line)
            if read_id is not None:
                yield read_id, "".join(seq_chunks)


def _compile_catalog_pattern(catalog):
    """Curated loci (e.g. FMR1) carry their own compound regex -- built the
    same way as strmie/scripts/pattern.py's htt_exact_match, i.e. from the
    locus's actual known repeat architecture, not a bare `(?:MOTIF)+`. See
    repeat_catalog.py for how each pattern was derived from real GRCh38
    reference sequence. Loci without one (e.g. C9orf72, which has no
    comparably-characterized interruption structure) fall back to a plain
    motif-repeat pattern."""
    if "pattern" in catalog:
        return re.compile(catalog["pattern"], re.I), catalog["repeat_unit_groups"]
    motif = catalog["motif"]
    return re.compile("(?:" + re.escape(motif.upper()) + ")+", re.I), None


def _units_in_match(match, repeat_unit_groups, motif_len):
    if repeat_unit_groups is None:
        return len(match.group(0)) // motif_len
    total = 0
    for group_name in repeat_unit_groups:
        span = match.group(group_name)
        if span:
            total += len(span) // motif_len
    return total


def longest_repeat_run(seq, catalog):
    """Best-scoring match of catalog's repeat pattern in seq, in repeat
    units. For curated loci this sums the repeat-unit groups across any
    interruption(s) (matching clinical convention, e.g. FMR1's total CGG
    count spanning its AGG interruption(s)); for plain motif loci it's the
    longest contiguous run. Returns None if the pattern doesn't match at
    all."""
    if not seq:
        return None
    pattern, repeat_unit_groups = _compile_catalog_pattern(catalog)
    motif_len = len(catalog["motif"])
    best = 0
    for m in pattern.finditer(seq.upper()):
        units = _units_in_match(m, repeat_unit_groups, motif_len)
        if units > best:
            best = units
    return best if best > 0 else None


def best_repeat_count(seq, catalog):
    """Best repeat-unit count checked in both the read's stored orientation
    and its reverse complement (unaligned reads have no reliable strand
    information), returning the larger of the two."""
    fwd = longest_repeat_run(seq, catalog)
    rev = longest_repeat_run(get_fast_rev_comp(seq), catalog)
    candidates = [x for x in (fwd, rev) if x is not None]
    return max(candidates) if candidates else None


def genotype_sample(fastq_path, locus, intorno=5, ii_threshold=False, ei_threshold=False):
    """Count `locus`'s repeat motif in every read of fastq_path and call
    the two most likely alleles. Returns a dict with the per-read counts
    DataFrame, the histogram, the (allele1, allele2) call (each may be
    "warning" if peak-calling failed, matching the HTT pipeline's own
    convention for an inconclusive sample), and the somatic Instability/
    Expansion indices (ii/ei; None if the alleles are a warning)."""
    catalog = load_catalog(locus) if isinstance(locus, str) else locus
    min_units = catalog.get("min_repeat_units", 2)

    def _warning(df):
        return {"catalog": catalog, "reads": df, "allele1": "warning", "allele2": "warning", "ii": None, "ei": None}

    records = []
    for read_id, seq in _read_fastq_or_fasta(fastq_path):
        count = best_repeat_count(seq, catalog)
        if count is not None:
            records.append((read_id, count))

    df = pd.DataFrame(records, columns=["ID", "Repeat_count"])
    if df.empty:
        return _warning(df)

    counts = df["Repeat_count"].value_counts()
    counts = counts[counts.index >= min_units]
    if counts.empty:
        return _warning(df)

    result = _two_highest_peaks(counts, intorno)
    if result is None:
        return _warning(df)

    allele1, allele2 = sorted(result)
    ii, ei = somatic_indices(df, allele1, allele2, ii_threshold=ii_threshold, ei_threshold=ei_threshold)
    return {"catalog": catalog, "reads": df, "allele1": allele1, "allele2": allele2, "ii": ii, "ei": ei}
