import argparse
import os

MOTIFS = {
    "canonical": "CAACAGCCGCCA",
    "m2": "CAGCAGCCGCCG",
    "loi_caa": "CAGCAGCCGCCA",
    "loi_cca": "CAACAGCCGCCG",
    "doi": "CAACAGCAACAGCCGCCA",
    "m6": "CAACAGCAACAGCCGCCG",
}
MOTIF_ORDER = ["canonical", "m2", "loi_caa", "loi_cca", "doi", "m6"]

SPLICE_START = 2875
SPLICE_END = 2969  # end of CCG repeat + "CCT" terminal in the native reference

# (sample_name, CAG1, CAG2, CCG1, CCG2, percentages [canonical, m2, loi_caa, loi_cca, doi, m6])
SAMPLES = [
    ("healthy_15_25_NOLOI",             15, 25, 7, 10, [100, 0, 0, 0, 0, 0]),
    ("healthy_17_22_NOLOI",             17, 22, 7, 10, [100, 0, 0, 0, 0, 0]),
    ("intermediate_12_34_NOLOI",        12, 34, 7, 10, [100, 0, 0, 0, 0, 0]),
    ("intermediate_16_31_DOI",          16, 31, 7, 10, [20, 0, 0, 0, 80, 0]),
    ("reduced_penetrance_17_37_LOICCA", 17, 37, 7, 10, [25, 0, 0, 75, 0, 0]),
    ("reduced_penetrance_19_39_LOI",    19, 39, 7, 10, [27, 0, 73, 0, 0, 0]),
    ("full_penetrance_16_39_DOI",       16, 39, 7, 10, [22, 0, 0, 0, 78, 0]),
    ("full_penetrance_23_48_LOI",       23, 48, 7, 10, [28, 0, 72, 0, 0, 0]),
    ("high_expansion_18_75_NOLOI",      18, 75, 7, 10, [100, 0, 0, 0, 0, 0]),
    ("high_expansion_18_84_NOLOI",      18, 84, 7, 10, [100, 0, 0, 0, 0, 0]),
    ("very_high_expansion_18_102_LOI",  18, 102, 7, 10, [29, 0, 71, 0, 0, 0]),
]


def load_ref(ref_path):
    seq = ""
    with open(ref_path) as f:
        for line in f:
            if not line.startswith(">"):
                seq += line.strip()
    return seq


def build_allele(seq, cag_len, ccg_len, motif_key):
    motif = MOTIFS[motif_key]
    repeat_region = ("CAG" * cag_len) + motif + ("CCG" * ccg_len) + "CCT"
    return seq[:SPLICE_START] + repeat_region + seq[SPLICE_END:]


def main():
    parser = argparse.ArgumentParser(
        description="Build engineered HTT exon 1 reference FASTAs for each Dataset 3 "
        "sample/interruption-motif combination, by splicing designed CAG/CCG/LOI-DOI "
        "repeat tracts into the native HTT locus sequence (chr4:3072000-3078000, hg38)."
    )
    parser.add_argument(
        "--ref",
        default=os.path.join(os.path.dirname(__file__), "htt_locus_chr4_3072000_3078000_hg38.fasta"),
    )
    parser.add_argument("--out-dir", default=os.path.join(os.path.dirname(__file__), "allele_refs"))
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    ref = load_ref(args.ref)
    manifest = []

    for name, cag1, cag2, ccg1, ccg2, pcts in SAMPLES:
        active_motifs = [(MOTIF_ORDER[i], pcts[i]) for i in range(6) if pcts[i] > 0]
        for motif_key, pct in active_motifs:
            allele1 = build_allele(ref, cag1, ccg1, motif_key)
            allele2 = build_allele(ref, cag2, ccg2, motif_key)
            fasta_path = os.path.join(args.out_dir, f"{name}__{motif_key}.fasta")
            with open(fasta_path, "w") as f:
                f.write(f">{name}_allele1_CAG{cag1}_CCG{ccg1}_{motif_key}\n{allele1}\n")
                f.write(f">{name}_allele2_CAG{cag2}_CCG{ccg2}_{motif_key}\n{allele2}\n")
            manifest.append((name, cag1, cag2, ccg1, ccg2, motif_key, pct, fasta_path))

    with open(os.path.join(args.out_dir, "manifest.tsv"), "w") as f:
        f.write("sample\tCAG1\tCAG2\tCCG1\tCCG2\tmotif\tpercent\tfasta_path\n")
        for row in manifest:
            f.write("\t".join(str(x) for x in row) + "\n")

    print(f"Built {len(manifest)} reference FASTA files for {len(SAMPLES)} samples -> {args.out_dir}")


if __name__ == "__main__":
    main()
