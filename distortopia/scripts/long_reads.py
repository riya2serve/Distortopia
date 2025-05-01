from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import re

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

    # ========== Visualization ==========
    print(f"Generating read visualization for: {out_path}")
    read_data = []
    for record in records:
        match = re.match(r"rep(\d+)_(.+?)_read(\d+)", record.id)
        if match:
            rep_id = int(match.group(1))
            contig = match.group(2)
            read_index = int(match.group(3))
            length = len(record.seq)
            read_data.append((rep_id, contig, read_index, length))

    df = pd.DataFrame(read_data, columns=["rep", "contig", "read_index", "length"])
    df = df.sort_values(by=["contig", "read_index"])

    fig, ax = plt.subplots(figsize=(12, 6))
    y_offset = 0
    yticks = []
    yticklabels = []

    for contig in df["contig"].unique():
        subset = df[df["contig"] == contig]
        for _, row in subset.iterrows():
            ax.barh(y=y_offset, width=row["length"], left=0, height=0.4, color="#4C72B0")
            y_offset += 1
        yticks.append(y_offset - len(subset) // 2)
        yticklabels.append(contig)
        y_offset += 1

    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels)
    ax.set_xlabel("Read Length (bp)")
    ax.set_title(f"Simulated Long Reads per Contig (Rep {rep})")

    os.makedirs("results", exist_ok=True)
    plot_path = os.path.join("results", f"long_read_map_rep{rep}.png")
    plt.tight_layout()
    plt.savefig(plot_path, dpi=300)
    plt.close()
    print(f"Visualization saved to: {plot_path}")


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



#=======
# EXAMPLE CLI
#=======
'''
Since .fasta file generated are huge, use the following function to check it worked
'''
#(base) riyarampalli@MacBookPro distortopia % head -n 10 genomes/f1_simulations/f1_reads_rep1.fasta

'''
Want to know total number of long reads simulated? Try this:
'''
#(base) riyarampalli@MacBookPro distortopia % grep ">" genomes/f1_simulations/f1_reads_rep1.fasta | wc -l