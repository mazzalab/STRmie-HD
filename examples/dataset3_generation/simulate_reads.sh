#!/bin/bash
# Simulate ONT long reads for each engineered Dataset 3 reference with Badread
# (v0.4.2), using the exact parameters reported in the manuscript Methods
# ("Dataset 3 synthetic generation"): R10.4.1 chemistry, mean read length
# 6000 bp (sd 1500), mean identity 98% (max 99.5%, sd 2), nanopore2023
# error/quality models, ~10,000 reads/sample, no junk/random/chimeric reads.
#
# Requires badread (pip install badread==0.4.2, or a dedicated conda env).
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MANIFEST="${1:-$SCRIPT_DIR/allele_refs/manifest.tsv}"
OUTDIR="${2:-$SCRIPT_DIR/dataset3_reads}"
WORKDIR="$SCRIPT_DIR/badread_raw"
mkdir -p "$OUTDIR" "$WORKDIR"

BP_PER_READ=4319
TARGET_READS=10000

# per-motif seed offsets, for deterministic regeneration
declare -A MOTIF_SEED=( [canonical]=1 [m2]=2 [loi_caa]=3 [loi_cca]=4 [doi]=5 [m6]=6 )

idx=0
tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r sample cag1 cag2 ccg1 ccg2 motif pct fasta; do
    idx=$((idx+1))
    quantity_bp=$(( BP_PER_READ * TARGET_READS * pct / 100 ))
    if [ "$quantity_bp" -lt 100000 ]; then
        quantity_bp=100000
    fi
    seed=$(( idx * 100 + ${MOTIF_SEED[$motif]} ))
    outfq="$WORKDIR/${sample}__${motif}.fastq"
    echo "=== $sample / $motif (pct=$pct, quantity=${quantity_bp}) seed=$seed ==="
    badread simulate --reference "$fasta" --quantity "${quantity_bp}" \
        --length 6000,1500 --identity 98,99.5,2 \
        --error_model nanopore2023 --qscore_model nanopore2023 \
        --junk_reads 0 --random_reads 0 --chimeras 0 --seed "$seed" \
        > "$outfq" 2> "$WORKDIR/${sample}__${motif}.log"
    n=$(( $(wc -l < "$outfq") / 4 ))
    echo "DONE $sample/$motif reads=$n"
done

echo "=== Concatenating per-sample FASTQs ==="
tail -n +2 "$MANIFEST" | awk -F'\t' '{print $1}' | sort -u | while read -r sample; do
    out="$OUTDIR/${sample}_1.fastq"
    rm -f "$out"
    for f in "$WORKDIR/${sample}__"*.fastq; do
        cat "$f" >> "$out"
    done
    gzip -f "$out"
    n=$(( $(zcat "${out}.gz" | wc -l) / 4 ))
    echo "SAMPLE $sample total_reads=$n"
done

echo "ALL_DATASET3_GENERATION_DONE"
