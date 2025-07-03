import random
import argparse
import numpy as np

def generate_repeat_sequence(cag_count, ccg_count, loi_pattern):
    """Generate a sequence with given CAG and CCG counts, and LOI pattern."""
    cag_part = "CAG" * cag_count
    ccg_part = "CCG" * ccg_count

    loi_dict = {
        "no": "CAACAGCCGCCA",
        "caacca": "CAGCAGCCGCCG",
        "caa": "CAGCAGCCGCCA",
        "cca": "CAACAGCCGCCG",
        "doi": "CAACAGCAACAGCCGCCA",
        "doicca": "CAACAGCAACAGCCGCCG"
    }

    loi_seq = loi_dict.get(loi_pattern, "CAACAGCCGCCA")

    return f"{cag_part}{loi_seq}{ccg_part}CCT"

def generate_fastq(output_file, cag1, cag2, ccg1, ccg2, loi_percentages, std_dev=2):
    """Generate a synthetic FASTQ file with reads."""
    total_reads = 10000  # Number of reads to generate
    read_length = 300

    # Generate CAG counts using a normal distribution centered on CAG1 and CAG2
    cag1_reads = np.random.normal(loc=cag1, scale=std_dev, size=int(total_reads * 0.4)).astype(int)
    cag2_reads = np.random.normal(loc=cag2, scale=std_dev, size=int(total_reads * 0.4)).astype(int)
    other_reads = np.random.randint(10, 100, size=int(total_reads * 0.2))

    # Combine all CAG values and clip to avoid negative values
    cag_values = np.concatenate([cag1_reads, cag2_reads, other_reads])
    cag_values = np.clip(cag_values, 1, 100)  # Ensure CAG values are within range

    random.shuffle(cag_values)

    # Generate exact number of reads for each LOI pattern
    loi_patterns = ["no", "caacca", "caa", "cca", "doi", "doicca"]
    loi_reads = []

    for pattern, percentage in zip(loi_patterns, loi_percentages):
        count = int(total_reads * (percentage / 100))
        for _ in range(count):
            cag_count = random.choice(cag_values)
            ccg_count = random.choice([ccg1, ccg2])
            seq = generate_repeat_sequence(cag_count, ccg_count, pattern)
            #seq = (seq[:read_length]).ljust(read_length, "A")
            loi_reads.append(f"@read{_}\n{seq}\n+\n{'I' * read_length}\n")

    # Shuffle the final reads to simulate randomness
    random.shuffle(loi_reads)

    # Write to FASTQ file
    with open(output_file, "w") as fq:
        fq.writelines(loi_reads)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--CAG1", type=int, default=17)
    parser.add_argument("--CAG2", type=int, default=49)
    parser.add_argument("--CCG1", type=int, default=7)
    parser.add_argument("--CCG2", type=int, default=10) # no, caacca, caa, cca, doi, doicca : ATTENZIONE ALLE CATEGORIE % NEL REPORT 
    parser.add_argument("--loi", nargs=6, type=int, default=[30, 20, 15, 10, 20, 5], help="Percentages for LOI patterns: no, caacca, caa, cca, doi, doicca")
    parser.add_argument("--output", type=str, default="synthetic_reads.fastq")
    parser.add_argument("--std_dev", type=float, default=2.0, help="Standard deviation for CAG distribution")

    args = parser.parse_args()

    generate_fastq(args.output, args.CAG1, args.CAG2, args.CCG1, args.CCG2, args.loi, args.std_dev)
    print(f"FASTQ file generated: {args.output}")

if __name__ == "__main__":
    main()
