#!/usr/bin/env python
# coding: utf-8
"""
Locus-aware BAM/CRAM read extraction for the --bam/--cram input modes.

Rather than teaching the (alignment-free, regex-based) calling logic to read
BAM/CRAM directly, this module extracts the reads relevant to the HTT CAG/CCG
locus into FASTQ, then hands off to the existing, unmodified FASTQ pipeline.

Extraction is index-based (fast: touches only the reads near the locus, not
the whole file) but is not a naive single-region query. A read supporting a
large repeat expansion is exactly the kind of read that tends to align
poorly (heavy soft-clipping, or the read/its mate marked unmapped), so a
plain "reads whose alignment position falls in the window" query would
silently lose the reads that matter most. To avoid that: any read in the
locus window whose mate maps outside the window is "rescued" by fetching
that exact mate position directly (via its recorded RNEXT/PNEXT, which BAM
indexing supports as a random-access lookup) -- the same anchoring principle
ExpansionHunter uses for in-repeat reads.

That still isn't enough on its own: a read entirely composed of repeat
sequence, with no unique flanking bases for the aligner to anchor on, can
end up with no alignment position at all -- and if its mate is unmapped
too (observed directly on real data: tens of thousands of such reads per
sample, all with an equally unmapped mate), there is no coordinate to
rescue by. These reads have no position-based signal to find them by, only
content. They are fetched via the indexed "no coordinate" bucket
(equivalent to `samtools view file.bam '*'`, which htslib still serves
without a full linear scan) and kept only if they match the locus by
content -- reusing STRmie-HD's own regex matcher, checked against both the
read's stored orientation and its reverse complement (its true orientation
is unknown without an alignment to derive it from).

The actual BAM -> FASTQ conversion is delegated to `samtools fastq`, which
already reverse-complements reverse-strand reads back to original sequencing
orientation. This matters because the standard (non --nanopore) calling path
has no reverse-complement fallback of its own -- it expects reads in a
single, consistent orientation.
"""

import os
import sys
import subprocess
import shutil

try:
    import pysam
except ImportError:  # pragma: no cover - surfaced as a clear runtime error, not at import time
    pysam = None

from strmie.scripts.pattern import htt_exact_match, get_fast_rev_comp


def _locus_match_needs_revcomp(seq):
    """Check whether seq matches STRmie-HD's own HTT CAG/CCG regex, in either
    orientation. Returns None if neither orientation matches (don't rescue
    this read), False if the stored orientation already matches, True if
    only the reverse complement does.

    This distinction matters because these are unmapped reads with no
    meaningful alignment, so their SAM "reverse strand" flag is essentially
    arbitrary -- but samtools fastq uses exactly that flag to decide whether
    to emit the stored sequence as-is or reverse-complemented. Getting a
    rescued read into the subset BAM isn't enough; its flag has to be set so
    the FASTQ it's converted to actually comes out in the orientation that
    matches, or the downstream regex-based caller will silently drop it too.
    """
    if not seq:
        return None
    if htt_exact_match(seq)[0] is not None:
        return False
    if htt_exact_match(get_fast_rev_comp(seq))[0] is not None:
        return True
    return None


# chr4 contig length is enough to disambiguate GRCh38 vs GRCh37/hg19 from a
# BAM/CRAM header alone, without requiring the user to state the build.
_BUILD_CHR4_LENGTH = {
    "grch38": 190214555,
    "grch37": 191154276,
}

# HTT exon 1 CAG+CCG repeat block coordinates, already established elsewhere
# in this project (ExpansionHunter variant catalogs, TRGT bed file).
_LOCUS_COORDS = {
    "grch38": {"repeat_start": 3074876, "repeat_end": 3074966},
    "grch37": {"repeat_start": 3076603, "repeat_end": 3076693},
}

# Flank sizes reproducing the 6 kb extraction window already validated
# elsewhere in this project (benchmarks/dataset5_v2/reference/htt_locus_repeat.bed,
# a locus reference spanning chr4:3072000-3078000 around this same repeat block).
_UPSTREAM_FLANK = 2876
_DOWNSTREAM_FLANK = 3034


def _resolve_samtools(samtools_bin):
    """Find the samtools executable. Honors an explicit absolute path, then
    PATH, then falls back to the bin/ directory of the currently running
    Python's own environment -- this keeps --bam/--cram self-contained even
    when the strmie entry point was invoked directly rather than via an
    activated conda environment (subprocess.run doesn't inherit an env's
    bin/ on PATH just because that env's python is what's executing)."""
    if os.path.isabs(samtools_bin) and os.path.isfile(samtools_bin):
        return samtools_bin
    found = shutil.which(samtools_bin)
    if found:
        return found
    candidate = os.path.join(sys.prefix, "bin", "samtools")
    if os.path.isfile(candidate):
        return candidate
    raise RuntimeError(
        f"The --bam/--cram input modes require 'samtools' (looked for "
        f"'{samtools_bin}' on PATH and in {os.path.join(sys.prefix, 'bin')}). "
        "Install it with `conda install -c bioconda samtools`."
    )


def _require_pysam():
    if pysam is None:
        raise ImportError(
            "The --bam/--cram input modes require the 'pysam' package, which is not "
            "installed in this environment. Install it with `conda install -c bioconda "
            "pysam` or `pip install pysam`, or convert your file to FASTQ yourself and "
            "use -f/--input instead."
        )


def _chrom_name(bam, prefer="chr4"):
    refs = set(bam.references)
    if prefer in refs:
        return prefer
    if "4" in refs:
        return "4"
    raise ValueError(
        "Could not find a 'chr4' or '4' contig in the BAM/CRAM header. "
        "Use --bam-region to specify the locus region manually (chrom:start-end)."
    )


def detect_build(bam):
    """Return 'grch38' or 'grch37' by matching the chr4/4 contig length in the
    header against known reference-build lengths. Raises ValueError if no
    match is found (caller should fall back to --locus-build/--bam-region)."""
    chrom = _chrom_name(bam)
    length = bam.lengths[bam.references.index(chrom)]
    for build, chrlen in _BUILD_CHR4_LENGTH.items():
        if length == chrlen:
            return build, chrom
    raise ValueError(
        f"Contig '{chrom}' has length {length}, which doesn't match a known "
        "GRCh38 (190214555) or GRCh37 (191154276) chr4 length. Auto-detection "
        "of the reference build failed -- pass --locus-build explicitly, or "
        "--bam-region chrom:start-end to bypass detection entirely."
    )


# A BAM/CRAM aligned against a reference that was itself already extracted
# down to just the repeat locus (a convention already used elsewhere in this
# project for ONT data, e.g. a 6 kb single-contig reference) needs no further
# windowing -- the whole contig already *is* the window. Treat any single,
# small contig this way rather than requiring it be named "chr4"/"4".
_SINGLE_CONTIG_LOCUS_MAX_LEN = 200_000


def resolve_region(bam, locus_build=None, region_override=None):
    """Return (chrom, start, end) [0-based, half-open, pysam convention]."""
    if region_override:
        return region_override

    if len(bam.references) == 1 and bam.lengths[0] <= _SINGLE_CONTIG_LOCUS_MAX_LEN:
        return bam.references[0], 0, bam.lengths[0]

    if locus_build and locus_build != "auto":
        chrom = _chrom_name(bam)
        build = locus_build
    else:
        build, chrom = detect_build(bam)

    c = _LOCUS_COORDS[build]
    start = max(0, c["repeat_start"] - _UPSTREAM_FLANK - 1)  # -1: 1-based -> 0-based
    end = c["repeat_end"] + _DOWNSTREAM_FLANK
    return chrom, start, end


def _mate_rescue_targets(reads, chrom, start, end):
    """Distinct (chrom, pos) pairs for mates of in-window reads that map
    outside the window and so would be missed by the window fetch alone."""
    targets = set()
    for r in reads:
        if not r.is_paired or r.mate_is_unmapped:
            continue
        mate_chrom = r.next_reference_name
        mate_pos = r.next_reference_start
        if mate_chrom is None or mate_pos is None:
            continue
        if mate_chrom != chrom or not (start <= mate_pos <= end):
            targets.add((mate_chrom, mate_pos))
    return targets


def extract_sample_to_fastq(
    path,
    sample_name,
    out_dir,
    reference=None,
    locus_build=None,
    region_override=None,
    samtools_bin="samtools",
    keep_subset_bam=False,
):
    """Extract HTT-locus-relevant reads from an indexed BAM/CRAM and write
    them to FASTQ(s) under out_dir. Returns {"paired": [r1, r2] or None,
    "single": path or None} -- both keys can be non-None at once (a sample
    can have both a proper R1/R2 pair and separately-rescued single-end
    content that needs to be folded into the same sample downstream, not
    discarded). Raises ValueError/ImportError/RuntimeError on failure.
    """
    _require_pysam()
    samtools_bin = _resolve_samtools(samtools_bin)

    is_cram = path.endswith(".cram")
    if is_cram and not reference:
        raise ValueError(f"--reference is required to decode CRAM file: {path}")

    os.makedirs(out_dir, exist_ok=True)

    mode = "rc" if is_cram else "rb"
    open_kwargs = {"reference_filename": reference} if reference else {}

    # pysam needs an index for random-access fetch(); build one next to the
    # input if missing and the directory is writable, matching how samtools
    # itself would require -X/an index for -L region queries.
    index_path = path + (".crai" if is_cram else ".bai")
    if not os.path.isfile(index_path):
        try:
            pysam.index(path)
        except Exception as exc:
            raise RuntimeError(
                f"{path} has no index ({index_path} not found) and could not be "
                f"auto-indexed ({exc}). Index it first with `samtools index {path}`."
            )

    with pysam.AlignmentFile(path, mode, **open_kwargs) as bam:
        chrom, start, end = resolve_region(bam, locus_build, region_override)
        window_reads = list(bam.fetch(chrom, start, end))
        rescue_targets = _mate_rescue_targets(window_reads, chrom, start, end)

        rescued_reads = []
        for mate_chrom, mate_pos in rescue_targets:
            rescued_reads.extend(bam.fetch(mate_chrom, max(0, mate_pos - 1), mate_pos + 1))

        # Reads with no alignment position at all (both mates unmapped, or
        # single-end with no mate to anchor on) can't be found by any
        # coordinate-based query. htslib still serves this "no coordinate"
        # bucket via the index (no full linear scan), same as `samtools view
        # file.bam '*'`. Keep only the ones that actually look like HTT reads
        # by content -- most of this bucket is unrelated background/off-target
        # sequence, not locus signal.
        content_rescued = []
        for r in bam.fetch("*"):
            if r.is_secondary or r.is_supplementary:
                continue
            needs_revcomp = _locus_match_needs_revcomp(r.query_sequence)
            if needs_revcomp is None:
                continue
            r.is_reverse = needs_revcomp
            content_rescued.append(r)

        subset_bam_path = os.path.join(out_dir, sample_name + ".htt_locus_subset.bam")
        seen = set()
        n_written = 0
        with pysam.AlignmentFile(subset_bam_path, "wb", template=bam) as out:
            for r in window_reads + rescued_reads + content_rescued:
                if r.is_secondary or r.is_supplementary:
                    continue
                key = (r.query_name, r.is_read1, r.is_read2)
                if key in seen:
                    continue
                seen.add(key)
                out.write(r)
                n_written += 1

    if n_written == 0:
        os.remove(subset_bam_path)
        raise RuntimeError(
            f"No reads found near the HTT locus ({chrom}:{start+1}-{end}) in {path}. "
            "If this file is aligned to a non-standard reference or a different "
            "assembly, pass --bam-region to specify the correct locus manually."
        )

    r1 = os.path.join(out_dir, sample_name + "_1.fastq.gz")
    r2 = os.path.join(out_dir, sample_name + "_2.fastq.gz")
    single = os.path.join(out_dir, sample_name + ".fastq.gz")
    # -0 (neither read1 nor read2 flag set -- true single-end reads) and -s
    # (a paired file's orphaned singletons) must NOT share a path: samtools
    # fastq opens each as its own gzip output stream, and two concurrent
    # writers to the same path interleave/corrupt the gzip container. Write
    # them separately and, if both end up with content, concatenate (valid
    # for gzip: concatenated gzip members decompress as one stream).
    unpaired = os.path.join(out_dir, sample_name + ".unpaired.fastq.gz")
    orphans = os.path.join(out_dir, sample_name + ".orphans.fastq.gz")

    cmd = [
        samtools_bin, "fastq",
        "-1", r1, "-2", r2,
        "-0", unpaired, "-s", orphans,
        subset_bam_path,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"samtools fastq failed for {path}:\n{result.stdout}\n{result.stderr}"
        )

    if not keep_subset_bam:
        os.remove(subset_bam_path)

    def _nonempty(p):
        # Checking the compressed file's raw byte size isn't reliable: an
        # empty gzip stream is itself ~35-40 bytes of header/footer, which
        # can sit right at a naive small-size threshold. Peek at whether
        # there's any decompressed content instead.
        if not os.path.isfile(p):
            return False
        import gzip as _gzip
        with _gzip.open(p, "rb") as fh:
            return fh.read(1) != b""

    have_unpaired = _nonempty(unpaired)
    have_orphans = _nonempty(orphans)
    if have_unpaired and have_orphans:
        with open(single, "wb") as out_f:
            for p in (unpaired, orphans):
                with open(p, "rb") as in_f:
                    shutil.copyfileobj(in_f, out_f)
        os.remove(unpaired)
        os.remove(orphans)
    elif have_unpaired:
        os.rename(unpaired, single)
        os.remove(orphans)
    elif have_orphans:
        os.rename(orphans, single)
        os.remove(unpaired)
    else:
        for p in (unpaired, orphans):
            if os.path.isfile(p):
                os.remove(p)

    # NOTE: content-rescued reads (see above) are matched by their own
    # sequence alone, with no requirement that their mate also matches --
    # most of them end up here as orphans/unpaired rather than in a proper
    # R1/R2 pair, since the mate typically covers unrelated flanking
    # sequence. This is commonly the *majority* of what rescue finds, not a
    # "handful of singletons" to discard alongside a paired sample: both
    # pieces are returned, and the caller is responsible for folding the
    # single-end content into the same sample rather than dropping it.
    paired = [r1, r2] if (_nonempty(r1) and _nonempty(r2)) else None
    if paired is None:
        for p in (r1, r2):
            if os.path.isfile(p):
                os.remove(p)
    single_out = single if (os.path.isfile(single) and _nonempty(single)) else None
    if single_out is None and os.path.isfile(single):
        os.remove(single)

    if paired is None and single_out is None:
        raise RuntimeError(f"samtools fastq produced no usable reads for {path}.")

    return {"paired": paired, "single": single_out}


def resolve_bam_cram_paths(input_paths, ext):
    """Mirror of utility.resolve_input_paths for .bam/.cram inputs: each
    entry may be a directory (all matching files inside are included) or a
    single file."""
    resolved = []
    for p in input_paths:
        if os.path.isdir(p):
            for name in sorted(os.listdir(p)):
                if name.endswith(ext) and os.path.isfile(os.path.join(p, name)):
                    resolved.append(os.path.join(p, name))
        elif os.path.isfile(p):
            if not p.endswith(ext):
                raise TypeError(f"Unsupported file format for input file: {p}. Expected {ext}")
            resolved.append(p)
        else:
            raise FileNotFoundError(f"Input path not found: {p}")

    if len(resolved) == 0:
        raise TypeError(f"No {ext} files found in the specified input path(s).")
    return resolved
