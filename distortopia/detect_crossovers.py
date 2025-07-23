import pandas as pd
import pysam
from collections import defaultdict

def load_snp_markers(marker_table_path):
    df = pd.read_csv(marker_table_path, sep="\t")

    # 🔥 Normalize CHROM format to match RefSeq (e.g., NC_003070.9)
    df["CHROM"] = df["CHROM"].astype(str)

    grouped = df.groupby("CHROM")
    snp_info = {
        chrom: set(positions_df["POS"].values)
        for chrom, positions_df in grouped
    }
    print("✅ SNP marker table loaded and grouped by chromosome.")
    return snp_info

def classify_read(read, snp_info):
    """Classify read by comparing its sequence to SNP alleles and detect crossover blocks."""
    read_seq = read.query_sequence
    read_start = read.reference_start
    read_chrom = read.reference_name
    read_pos = read.get_reference_positions()

    markers = snp_info.get(read_chrom)
    if not markers:
        return None

    pattern = []
    crossover_positions = []
    prev = None

    for snp_pos, ref_base, alt_base in zip(markers["positions"], markers["ref_alleles"], markers["alt_alleles"]):
        if snp_pos in read_pos:
            read_index = read_pos.index(snp_pos)
            base = read_seq[read_index]
            if base == ref_base:
                allele = "T"
            elif base == alt_base:
                allele = "L"
            else:
                allele = "N"
            pattern.append(allele)

            if prev and allele != prev and allele in "TL" and prev in "TL":
                crossover_positions.append((read.reference_name, snp_pos))
            prev = allele

    if not pattern:
        return None

    print(f"🔬 Pattern for read {read.query_name}: {''.join(pattern)} with {len(crossover_positions)} crossovers")

    return {
        "qname": read.query_name,
        "chrom": read.reference_name,
        "start": read.reference_start,
        "end": read.reference_end,
        "length": read.reference_length,
        "pattern": "".join(pattern),
        "crossovers": crossover_positions,
    }

def detect_crossovers(bam_path, marker_table_path, output_path="crossovers.tsv"):
    snp_info = load_snp_markers(marker_table_path)
    bam = pysam.AlignmentFile(bam_path, "rb")

    results = []
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        classification = classify_read(read, snp_info)
        if classification:
            results.append(classification)

    bam.close()

    print(f"📦 Total classified reads: {len(results)}")

    df = pd.DataFrame(results)

    if not df.empty and "crossovers" in df.columns:
        df["num_crossovers"] = df["crossovers"].apply(lambda x: len(x))
    else:
        print("⚠️ No reads with crossover information found.")
        df["num_crossovers"] = 0

    df.to_csv(output_path, sep="\t", index=False)
    print(f"✅ Crossover detection complete. Results saved to {output_path}.")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Detect crossovers from BAM file using SNP marker table.")
    parser.add_argument("--bam", required=True, help="Path to F1 BAM file.")
    parser.add_argument("--markers", required=True, help="Path to SNP marker table (from compare_variants.py).")
    parser.add_argument("--out", default="crossovers.tsv", help="Output path for crossover TSV.")
    args = parser.parse_args()

    detect_crossovers(args.bam, args.markers, args.out)

