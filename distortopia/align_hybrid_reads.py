import subprocess
import os
import sys

def log(message, log_box=None):
    print(message)
    if log_box:
        log_box.write(message)

def align_reads(fq_gz, parent_ref, output_prefix, log_box=None):
    bam = f"{output_prefix}.bam"
    sorted_bam = f"{output_prefix}.sort.bam"

    log(f"🔗 Aligning to {parent_ref}...", log_box)
    cmd_align = ["minimap2", "-ax", "map-pb", parent_ref, fq_gz]
    cmd_view = ["samtools", "view", "-Sb", "-"]
    cmd_sort = ["samtools", "sort", "-", "-o", sorted_bam]
    cmd_index = ["samtools", "index", sorted_bam]

    try:
        # Run minimap2 -> samtools view -> samtools sort
        align = subprocess.Popen(cmd_align, stdout=subprocess.PIPE)
        view = subprocess.Popen(cmd_view, stdin=align.stdout, stdout=subprocess.PIPE)
        sort = subprocess.Popen(cmd_sort, stdin=view.stdout)
        sort.communicate()

        # Index the sorted BAM
        subprocess.run(cmd_index, check=True)

        log(f"✅ Alignment complete: {sorted_bam}", log_box)
        return sorted_bam
    except Exception as e:
        log(f"❌ Alignment failed: {e}", log_box)
        raise RuntimeError("Alignment pipeline failed.") from e

# For CLI use (not Streamlit)
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--fq", required=True, help="Hybrid .fq.gz file")
    parser.add_argument("--ref1", required=True, help="Parent 1 reference")
    parser.add_argument("--ref2", required=True, help="Parent 2 reference")
    args = parser.parse_args()

    align_reads(args.fq, args.ref1, "F1_to_thaliana")
    align_reads(args.fq, args.ref2, "F1_to_lyrata")

