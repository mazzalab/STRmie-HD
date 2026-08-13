#!/usr/bin/env python
# coding: utf-8
"""
Somatic instability indices for the generic (non-HTT) submodule.

Reuses strmie.scripts.indices' instabilityIndex/expansionIndex exactly as
they are for the HTT pipeline -- imported, not copied or modified. Those
functions expect a distribution dataframe shaped like
strmie.scripts.utility.create_df_distribution's output: one row per observed
repeat length, with a "count" column (reads at that length) and a
"CAG_repeat" column (the length itself, despite the HTT-specific column
name). That shape is really just "a repeat-length histogram" -- nothing
about the computation itself is CAG-specific -- so building the same shape
from this submodule's generic Repeat_count column is legitimate, zero-risk
reuse: strmie/scripts/indices.py is not touched.

Deliberately NOT reused here:
  - histogramRatioIndex: keyed to a fixed clinical CAG cutpoint specific to
    HD, with no equivalent for FMR1/C9orf72.
  - LOI/LOI_CCA/DOI interruption-percentage calculations: HTT's specific
    interruption taxonomy. FMR1's own interruption signal (AGG loss) is
    instead folded directly into the curated repeat count itself (see
    repeat_catalog.py), not tracked as a separate metric here.
"""

from strmie.scripts.indices import instabilityIndex, expansionIndex


def _distribution_df(reads_df):
    dist = reads_df["Repeat_count"].value_counts().to_frame()
    dist["CAG_repeat"] = dist.index
    return dist


def somatic_indices(reads_df, allele1, allele2, ii_threshold=False, ei_threshold=False):
    """(instability_index, expansion_index) for this sample's called
    alleles, or (None, None) if the alleles are warnings/too close together
    to compute (mirrors indices.py's own convention of returning a message
    string in that case, converted here to None for JSON/Excel-friendliness)."""
    if allele1 in (None, "warning") or allele2 in (None, "warning"):
        return None, None
    dist = _distribution_df(reads_df)
    ii = instabilityIndex(dist, allele1, allele2, pcrFiltering=ii_threshold)
    ei = expansionIndex(dist, allele1, allele2, pcrFiltering=ei_threshold)
    ii = ii if isinstance(ii, (int, float)) else None
    ei = ei if isinstance(ei, (int, float)) else None
    return ii, ei
