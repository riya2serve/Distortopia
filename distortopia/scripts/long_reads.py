rom Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import os
import argparse
import numpy as np

def load_f1_fasta(fasta_path):
    return {record.id: str(record.seq) for record in SeqIO.parse(fasta_path, "fasta")}

def simulate_long_reads(f1_genome, f1_table, rep, output_dir, read_len=15000, coverage=60):
    records = []
    subtable = f1_table[f1_table["rep"] == rep]
    seen_contigs = set(subtable["ref_chrom"].unique())
    
    for contig_id in seen_contigs:
        hybrid_id = f"{contig_id}_rep{rep}"
        if hybrid_id not in f1_genome:
            continue
        seq = f1_genome[hybrid_id]
        seq_len = len(seq)
        num_reads = (seq_len * coverage) // read_len
        
        for i in range(num_reads):
            start = np.random.randint(0, seq_len - read_len + 1)
            fragment = seq[start:start + read_len]
            read_id = f"rep{rep}_{contig_id}_read{i}"
            records.append(SeqRecord(Seq(fragment), id=read_id, description="simulated_long_read"))
    
    out_path = os.path.join(output_dir, f"f1_reads_rep{rep}.fasta")
    with open(out_path, "w") as out_handle:
        SeqIO.write(records, out_handle, "fasta-2line")

def main(ref_fasta, alt_fasta, snp_tsv, output_dir, n_reps):
    os.makedirs(output_dir, exist_ok=True)
    f1_table = pd.read_csv(os.path.join(output_dir, "f1_table.csv"))
    
    for rep in range(1, n_reps + 1):
        print(f"Simulating long reads for rep {rep}...")
        fasta_path = os.path.join(output_dir, f"f1_genome_rep{rep}.fna")
        f1_genome = load_f1_fasta(fasta_path)
        simulate_long_reads(f1_genome, f1_table, rep, output_dir)

def parse_args():
    parser = argparse.ArgumentParser(description="Simulate long reads from F1 hybrid genomes.")
    parser.add_argument("--ref-fasta", required=True, help="Path to parent 1 FASTA (not used but kept for compatibility)")
    parser.add_argument("--alt-fasta", required=True, help="Path to parent 2 FASTA (not used but kept for compatibility)")
    parser.add_argument("--snp-tsv", required=True, help="Path to SNP table (not used but kept for compatibility)")
    parser.add_argument("--outdir", default="genomes/f1_simulations", help="Directory where F1 genome files are stored")
    parser.add_argument("--reps", type=int, default=1, help="Number of F1 replicates")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    main(args.ref_fasta, args.alt_fasta, args.snp_tsv, args.outdir, args.reps)