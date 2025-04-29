from Bio import SeqIO #BioPython module for reading and writing sequence files 
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord
import pandas as pd #for working with dataframes and generated HTML style summaries
import os #interacts with operating system 
import argparse #allows for command-line arguments/flags from/for users

"""
This script loads A. thaliana and A. lyrata .fna files and then parses the SNP alignment summary.
It allows users to extract both reference and target bases at SNP positions. 
With this, users can simulate a haploid F1 hybrid, with target SNPs incorporated into the reference
genome. 
"""
def load_fasta(fasta_path):
    #loads query ('reference') and target FASTA sequences into dictionaries
    return {record.id: list(str(record.seq)) for record in SeqIO.parse(fasta_path, "fasta")}

def parse_snps(csv_path):
    df = pd.read_csv(csv_path, sep="\t")
    snp_data = []

    for _, row in df.iterrows():
        query_id = row["Query"].split()[0]  # e.g., BASP01000003.1
        target_id = row["Target"].split()[0]  # e.g., CP087130.2
        positions = row["SNP_Positions"]
        if pd.notna(positions):
            snps = [int(pos.strip()) for pos in str(positions).split(",")]
            snp_data.append((query_id, target_id, snps))
    
    return snp_data

def simulate_hybrid(ref_genome, tgt_genome, snp_data):
    #generating hybrid by replacing ref. base with target base
    hybrid = ref_genome.copy()

    for query_id, target_id, snps in snp_data:
        if query_id not in ref_genome or target_id not in tgt_genome:
            continue  # Skip contigs not found

        ref_seq = ref_genome[query_id]
        tgt_seq = tgt_genome[target_id]

        for pos in snps:
            idx = pos - 1  # Convert to 0-based
            try:
                hybrid[query_id][idx] = tgt_seq[idx]  # Replace with target base
            except IndexError:
                continue  # Handle edge case if index exceeds contig length

    return hybrid

def write_fasta(output_dict, output_path):
    #writes hybrid sequence to a FASTA file
    records = []
    for contig_id, seq in output_dict.items():
        record = SeqRecord(Seq("".join(seq)), id=contig_id, description="hybrid")
        records.append(record)
    SeqIO.write(records, output_path, "fasta")

def parse_args():
    """
    Command-line argument interface.
    """
    parser = argparse.ArgumentParser(description="Simulate pseudo-F1 hybrid FASTA from two parent genomes.")
    parser.add_argument("--ref-dir", required=True, help="Path to A. thaliana reference .fna file")
    parser.add_argument("--query-dir", required=True, help="Path to A. lyrata target .fna file")
    parser.add_argument("--snp", required=True, help="Path to SNP summary table (TSV)")
    parser.add_argument("--out", default="pseudo_F1_hybrid.fasta", help="Output FASTA filename")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()  # <--- call and assign it here

    print("Loading genomes...")
    ref_genome = load_fasta(args.ref_dir)
    tgt_genome = load_fasta(args.query_dir)

    print("Parsing SNP table...")
    snp_data = parse_snps(args.snp)

    print("Simulating F1 hybrid...")
    hybrid_genome = simulate_hybrid(ref_genome, tgt_genome, snp_data)

    print(f"Writing output to {args.out}")
    write_fasta(hybrid_genome, args.out)


## ==========
# EXAMPLE CLI 
## ==========

#bash
##python simf1poly.py \
  ###--ref-dir user-data/Arabidopsis_thaliana/ncbi_dataset/data/GCA_000001735.2/GCA_000001735.2_TAIR10.1_genomic.fna \
  ###--query-dir user-data/Arabidopsis_lyrata/ncbi_dataset/data/GCA_000004255.1/GCA_000004255.1_v.1.0_genomic.fna \
  ###--snp snp_positions.tsv \
  ###--out F1_hybrid.fna
