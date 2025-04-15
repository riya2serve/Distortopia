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

def compare_fastas(ref_fasta, query_fasta, output_vcf):
    """
    Compare two FASTA files base-by-base and output SNPs to VCF.
    """
    ref_seqs = SeqIO.to_dict(SeqIO.parse(ref_fasta, "fasta"))
    query_seqs = SeqIO.to_dict(SeqIO.parse(query_fasta, "fasta"))

    with open(output_vcf, "w") as vcf:
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write("##source=BiopythonSNPComparer\n")
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for chrom in ref_seqs:
            if chrom not in query_seqs:
                print(f"[SKIP] {chrom} not found in query genome.")
                continue

            ref_seq = str(ref_seqs[chrom].seq).upper()
            query_seq = str(query_seqs[chrom].seq).upper()

            min_len = min(len(ref_seq), len(query_seq))
            snp_count = 0

            for i in range(min_len):
                r = ref_seq[i]
                q = query_seq[i]
                if r != q and r in "ACGT" and q in "ACGT":
                    vcf.write(f"{chrom}\t{i+1}\t.\t{r}\t{q}\t.\tPASS\t.\n")
                    snp_count += 1

            print(f"Compared {chrom} → {snp_count:,} SNPs found.")

    print(f"\nSNP comparison complete → {output_vcf}")

def parse_args():
    parser = argparse.ArgumentParser(description="Interactively compare two genome FASTAs and identify SNPs")
    parser.add_argument("--ref-dir", required=True, help="Directory containing FASTA(s) for species 1")
    parser.add_argument("--query-dir", required=True, help="Directory containing FASTA(s) for species 2")
    parser.add_argument("--out", required=True, help="Output VCF path")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    # Interactive FASTA selection
    ref_fasta = choose_fasta("Species 1", args.ref_dir)
    query_fasta = choose_fasta("Species 2", args.query_dir)

    # Compare and output SNPs
    compare_fastas(ref_fasta, query_fasta, args.out)


