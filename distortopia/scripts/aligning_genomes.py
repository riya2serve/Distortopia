import subprocess
import argparse
import os

"""
This script allows users to align their F1 hybrid genome to a reference genome. 
It accepts two FASTA files and uses minimap2 to in .paf mode.
It outputs an alignment file (.paf)
"""

def run_minimap2(reference_fasta, hybrid_fasta, output_file, mode="paf"):
    if mode == "paf":
        output_format = "paf"
        minimap2_args = ["-x", "asm5"]
    elif mode == "sam":
        output_format = "sam"
        minimap2_args = ["-ax", "asm5"]
    else:
        raise ValueError("Mode must be 'paf' or 'sam'.")

    cmd = ["minimap2"] + minimap2_args + [reference_fasta, hybrid_fasta]

    with open(output_file, "w") as out_f:
        subprocess.run(cmd, stdout=out_f)

    print(f"Minimap2 alignment complete. Output saved to {output_file}")

def parse_args():
    parser = argparse.ArgumentParser(description="Align hybrid FASTA to reference using Minimap2.")
    parser.add_argument("--ref-dir", required=True, help="Reference genome FASTA (e.g., A. thaliana)")
    parser.add_argument("--query-dir", required=True, help="Query genome FASTA (e.g., F1 hybrid)")
    parser.add_argument("--out", required=True, help="Output alignment file (e.g., hybrid_vs_ref.paf or .sam)")
    parser.add_argument("--mode", choices=["paf", "sam"], default="paf", help="Output format: paf (default) or sam")
    return parser.parse_args()

if __name__ == "__main__":
	args = parse_args()  # <--- call and assign it here

	run_minimap2(args.ref_dir, args.query_dir, args.out, args.mode)

## ==========
# EXAMPLE CLI 
## ==========
"""
User is in PAF mode and is using it w/ argparse flags. 
Note: SAM mode is useful if user is seeking to view hybrid genome using IGV
"""
#bash
##python aligning_genomes.py \
  ###--ref user-data/Arabidopsis_thaliana/ncbi_dataset/data/GCA_000001735.2/GCA_000001735.2_TAIR10.1_genomic.fna \
  ###--query F1_hybrid.fna \
  ###--out hybrid_vs_ref.paf

"""
Users can run this script as many times, depending on how many reference genomes they want to compare
their F1 haploid genomes to. Generated .paf files should have different names, so as to not overwrite each time.

To visualize alignment, users should navigate to <https://dgenies.toulouse.inra.fr>. 
D-genies is a free, open-source tool to compare two genomes. It supports large genomes and allows
users to interact with the dot plot to improve the visualisation.
"""





