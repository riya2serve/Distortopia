import os
import glob
from Bio import SeqIO
import argparse

"""
This script will take two input FASTA files (--ref-dir and --query-dir).
It will align the query to the reference using 'minimap2.'
It will then call variants using paftools.js
Variant calls will be ouput onto a real VCF, which can be used later
to simulate F1 hybrid genotypes or detect segregation distortion."
"""

def choose_fasta(species_name, species_dir):
    """
    Search for FASTA files in a species directory and prompt user to select one.
    """
    fasta_paths = glob.glob(f"{species_dir}/**/*genomic.fna", recursive=True)
    if not fasta_paths:
        raise FileNotFoundError(f"No FASTA files found for {species_name} in {species_dir}")

    print(f"\n[SELECT FASTA FOR {species_name}]")
    for i, path in enumerate(fasta_paths):
        print(f"{i+1}. {path}")

    while True:
        try:
            index = int(input(f"Enter number [1-{len(fasta_paths)}]: ")) - 1
            if 0 <= index < len(fasta_paths):
                return fasta_paths[index]
            else:
                print("Invalid selection. Try again.")
        except ValueError:
            print("Please enter a valid number.")

def align_call_variants(ref_fasta, query_fasta, output_vcf, paftools_path):
    """
    Align query to reference using minimap2 and call variants with paftools.js.
    """
    paf_file = output_vcf.replace(".vcf", ".paf")

    print(f"\nAligning query → reference using minimap2...")
    with open(paf_file, "w") as paf_out:
        subprocess.run(["minimap2", "-x", "asm5", ref_fasta, query_fasta],
                       stdout=paf_out, check=True)

    print(f"Calling variants with paftools.js...")
    with open(output_vcf, "w") as vcf_out:
        subprocess.run(["node", paftools_path, "call", paf_file],
                       stdout=vcf_out, check=True)

    print(f"\nVariant calling complete → {output_vcf}")

def parse_args():
    parser = argparse.ArgumentParser(description="Align two genome FASTAs and call real variants")
    parser.add_argument("--ref-dir", required=True, help="Directory containing FASTA(s) for reference species")
    parser.add_argument("--query-dir", required=True, help="Directory containing FASTA(s) for query species")
    parser.add_argument("--out", required=True, help="Output VCF path")
    parser.add_argument("--paftools", default="minimap2/paftools.js", help="Path to paftools.js")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    # Interactive FASTA selection
    ref_fasta = choose_fasta("Species 1", args.ref_dir)
    query_fasta = choose_fasta("Species 2", args.query_dir)

    # Compare and output SNPs
    align_and_call_variants(ref_fasta, query_fasta, args.out, args.paftools)



