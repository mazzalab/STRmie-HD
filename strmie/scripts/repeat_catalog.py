#!/usr/bin/env python
# coding: utf-8
"""
Repeat catalog for the generic (non-HTT) submodule.

This is intentionally separate from the HTT-specific pipeline: it does not
import from, modify, or get imported by strmie/main.py or
strmie/scripts/pattern.py's HTT-specific functions. It demonstrates that
STRmie-HD's core alignment-free, regex-based counting approach generalizes
to other tandem repeat loci by defining them data-natively (motif +
reference region + optional known size ranges), rather than by writing new
per-disease code.

A catalog entry needs only a repeat motif to be usable for counting; the
reference region is included because it can be reused as-is by
strmie/scripts/bam_extract.py's --bam/--cram extraction (already
region/build-parametrized), and the size ranges are included purely as
reference context for interpreting results, not used by the matching logic
itself.

Two optional fields let a locus carry a curated, interruption-aware pattern
instead of a bare motif-repeat, the same way HTT's own htt_exact_match is
curated from HTT's real repeat architecture rather than being a generic
`(?:CAG)+`:
  - "pattern": a regex string (case-insensitive) with named groups for each
    contiguous repeat-unit block.
  - "repeat_unit_groups": the group names to sum (in repeat units) for the
    total reported count.
A locus with neither field falls back to a plain `(?:MOTIF)+` longest-run
match -- the honest choice when, unlike HTT or FMR1, the locus has no
comparably well-characterized interruption architecture in the literature.

Two further optional fields, "flank_upstream"/"flank_downstream", enable an
opt-in fuzzy-matching mode (see generic_repeat.py's fuzzy_repeat_count),
built on the same design idea as strmie/scripts/pattern.py's own
htt_nanopore_match: fuzzy-match short, non-repetitive sequences immediately
flanking the repeat tract (tolerant of a bounded edit distance, via the
`regex` module), then estimate the repeat-unit count from the length of the
region between them, rather than requiring every base inside the tract to
exactly spell the motif -- tolerant of read errors landing inside the
repeat itself, not just at noisy read termini. This is a fresh, separate
implementation for arbitrary catalog loci, not a reuse of
htt_nanopore_match itself (whose LOI/DOI interruption-motif-block search is
HTT-specific); strmie/scripts/pattern.py is not imported from or modified
here. Exact matching (pattern/motif above) remains the default for every
locus; fuzzy matching is opt-in per run via strmie-repeat's --nanopore flag
and only usable for a locus that defines these two flanks.
"""

import json
import os

# Built-in catalog entries, beyond HTT (which has its own dedicated,
# interruption-motif-aware pipeline in strmie/main.py and is intentionally
# not duplicated here). GRCh38 coordinates, 1-based inclusive.
BUILTIN_CATALOG = {
    "FMR1": {
        "gene": "FMR1",
        "disease": "Fragile X Syndrome",
        "motif": "CGG",
        # Curated, interruption-aware pattern -- built the same way as
        # strmie/scripts/pattern.py's htt_exact_match: derived from the
        # locus's real reference sequence (verified against GRCh38
        # chrX:147911950-147912200, which shows (CGG)10 AGG (CGG)9), not a
        # bare `(?:CGG)+`. FMR1's CGG array is normally interrupted by up to
        # two AGG triplets roughly every ~9-10 repeats; clinically reported
        # CGG size counts the whole array across interruption(s), and loss
        # of these AGG interruptions is itself the established marker of
        # instability/expansion risk for this locus -- directly analogous
        # to HTT's LOI-CAA concept, just embedded within the tract instead
        # of at its boundary.
        "pattern": r"(?P<CGG_1>(cgg)+)(?:agg(?P<CGG_2>(cgg)+))?(?:agg(?P<CGG_3>(cgg)+))?",
        "repeat_unit_groups": ["CGG_1", "CGG_2", "CGG_3"],
        # 20bp immediately outside the repeat tract, read directly from
        # GRCh38 chrX:147911950-147912200 (repeat span verified at
        # chrX:147912051-147912110, matching "region" below almost exactly).
        "flank_upstream": "CCAGGGGGCGTGCGGCAGCG",
        "flank_downstream": "CTGGGCCTCGAGCGCCCGCA",
        "min_repeat_units": 3,
        "reference_build": "GRCh38",
        "region": "chrX:147912050-147912110",
        "normal_range": [5, 44],
        "premutation_range": [55, 200],
        "full_mutation_range": [200, None],
    },
    "C9orf72": {
        "gene": "C9orf72",
        "disease": "ALS/FTD",
        # No "pattern" field: verified against real reference (GRCh38
        # chr9:27573400-27573650, minus strand -- C9orf72 is minus-strand;
        # the GGGGCC array reads cleanly on that strand: GGGGCC-GGGGCC-
        # GGGGCC-...), but unlike FMR1/HTT, this locus has no comparably
        # well-characterized internal interruption motif in the literature
        # to curate a compound pattern from -- pathogenic alleles are
        # essentially one long uninterrupted GGGGCC run, and clinical
        # sizing (repeat-primed PCR) treats it the same way. A plain
        # longest-run match is the accurate, not merely convenient, choice
        # here.
        "motif": "GGGGCC",
        # 20bp immediately outside the repeat tract, read directly from
        # GRCh38 chr9:27573400-27573650, gene-sense (minus/revcomp) strand,
        # same orientation as the "motif" above.
        "flank_upstream": "AACTCAGGAGTCGCGCGCTA",
        "flank_downstream": "GGGGCGTGGTCGGGGCGGGC",
        "min_repeat_units": 2,
        "reference_build": "GRCh38",
        "region": "chr9:27573528-27573546",
        "normal_range": [2, 23],
        "pathogenic_range": [30, None],
    },
}


_CATALOG_BY_UPPER_NAME = {name.upper(): name for name in BUILTIN_CATALOG}


def load_catalog(name_or_path):
    """Resolve a catalog entry: a built-in locus name (case-insensitive),
    or a path to a user-supplied JSON file with the same schema as the
    built-in entries above."""
    canonical = _CATALOG_BY_UPPER_NAME.get(name_or_path.upper())
    if canonical is not None:
        return dict(BUILTIN_CATALOG[canonical])

    if os.path.isfile(name_or_path):
        with open(name_or_path) as f:
            entry = json.load(f)
        required = {"gene", "motif"}
        missing = required - set(entry)
        if missing:
            raise ValueError(f"Catalog file {name_or_path} is missing required field(s): {missing}")
        entry.setdefault("min_repeat_units", 2)
        return entry

    raise ValueError(
        f"Unknown locus/catalog '{name_or_path}'. Built-in loci: {', '.join(BUILTIN_CATALOG)}. "
        "Otherwise, pass a path to a JSON catalog file with at least 'gene' and 'motif' fields."
    )


def list_builtin_loci():
    return list(BUILTIN_CATALOG)
