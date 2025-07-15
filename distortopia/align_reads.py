import subprocess
import os

def run_alignment(reads_path, reference_fasta, output_prefix, aligner='minimap2'):
    """
    Align long reads to a reference genome and generate a sorted, indexed BAM file.

    Args:
        reads_path (str): Path to simulated long reads (.fq or .fq.gz)
        reference_fasta (str): Path to reference genome FASTA
        output_prefix (str): Prefix for output files (e.g. "sim_thaliana")
        aligner (str): Aligner to use (currently supports 'minimap2')
    """
    # Step 1: Build index (.mmi)
    ref_index = f"{reference_fasta}.mmi"
    if not os.path.exists(ref_index):
        print(f"[INFO] Building minimap2 index for {reference_fasta}")
        subprocess.run([aligner, "-d", ref_index, reference_fasta], check=True)

    # Step 2: Align reads to reference
    sam_file = f"{output_prefix}.sam"
    print(f"[INFO] Aligning reads with minimap2: {reads_path} → {sam_file}")
    with open(sam_file, "w") as sam_out:
        subprocess.run([
            aligner, "-ax", "map-pb", ref_index, reads_path
        ], stdout=sam_out, check=True)

    # Step 3: Convert to BAM, sort, and index
    bam_file = f"{output_prefix}.bam"
    sorted_bam = f"{output_prefix}.sort.bam"

    print(f"[INFO] Converting SAM to BAM: {sam_file} → {bam_file}")
    subprocess.run(["samtools", "view", "-Sb", sam_file], stdout=open(bam_file, "wb"), check=True)

    print(f"[INFO] Sorting BAM: {bam_file} → {sorted_bam}")
    subprocess.run(["samtools", "sort", "-o", sorted_bam, bam_file], check=True)

    print(f"[INFO] Indexing BAM: {sorted_bam}")
    subprocess.run(["samtools", "index", sorted_bam], check=True)

    print(f"[DONE] Alignment complete for {output_prefix}.")


# Optional command-line usage
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Align simulated long reads to a reference genome.")
    parser.add_argument("--reads", required=True, help="Path to .fq or .fq.gz file")
    parser.add_argument("--ref", required=True, help="Path to reference .fna file")
    parser.add_argument("--prefix", required=True, help="Output file prefix (e.g., sim_thaliana)")

    args = parser.parse_args()
    run_alignment(args.reads, args.ref, args.prefix)

