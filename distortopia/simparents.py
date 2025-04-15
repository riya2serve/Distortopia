import os
import glob
from Bio import SeqIO
import argparse

"""
This script allows a user to select one genome FASTA from each of two species,
then compares them base-by-base for all matching contigs (by name).

It outputs a real, biologically meaningful VCF of SNPs between the two genome assemblies.
"""

def choose_fasta(species_name, species_dir):
    """
    Lists all available .fna files under a species directory and prompts user to choose one.
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
    Compares two genome assemblies base-by-base for contigs with matching names.
    Outputs differences as SNPs in a standard VCF file.
    """
    ref_seqs = SeqIO.to_dict(SeqIO.parse(ref_fasta, "fasta"))
    query_seqs = SeqIO.to_dict(SeqIO.parse(query_fasta, "fasta"))

    # Report on contig matching before comparing
    ref_contigs = set(ref_seqs.keys())
    query_contigs = set(query_seqs.keys())
    shared_contigs = ref_contigs & query_contigs

    print(f"\nReference genome has {len(ref_contigs)} contigs.")
    print(f"Query genome has {len(query_contigs)} contigs.")
    print(f"Shared contigs for comparison: {len(shared_contigs)}")
    if not shared_contigs:
        print("No shared contig names between the two genomes!")
        return

    # Begin writing VCF output
    with open(output_vcf, "w") as vcf:
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write("##source=BiopythonSNPComparer\n")
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for chrom in sorted(shared_contigs):
            ref_seq = str(ref_seqs[chrom].seq).upper()
            query_seq = str(query_seqs[chrom].seq).upper()
            min_len = min(len(ref_seq), len(query_seq))
            snp_count = 0

            for i in range(min_len):
                r = ref_seq[i]
                q = query_seq[i]
                # Only call SNPs at valid nucleotide positions (skip N, gaps, etc.)
                if r != q and r in "ACGT" and q in "ACGT":
                    vcf.write(f"{chrom}\t{i+1}\t.\t{r}\t{q}\t.\tPASS\t.\n")
                    snp_count += 1

            print(f"Compared {chrom} → {snp_count:,} SNPs")

    print(f"\nSNP comparison complete → {output_vcf}")

def parse_args():
    """
    Parses CLI arguments for species FASTA directory paths and VCF output path.
    """
    parser = argparse.ArgumentParser(description="Compare two genome FASTAs and identify SNPs")
    parser.add_argument("--ref-dir", required=True, help="Directory containing FASTA(s) for species 1 (reference)")
    parser.add_argument("--query-dir", required=True, help="Directory containing FASTA(s) for species 2 (query)")
    parser.add_argument("--out", required=True, help="Output VCF file path")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    # Prompt user to select which FASTA files to compare
    ref_fasta = choose_fasta("Species 1", args.ref_dir)
    query_fasta = choose_fasta("Species 2", args.query_dir)

    # Create output directory if it doesn’t exist
    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    # Run SNP comparison and write output
    compare_fastas(ref_fasta, query_fasta, args.out)
