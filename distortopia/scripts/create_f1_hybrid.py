from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import os
import argparse
import numpy as np

def load_fasta(fasta_path):
    """
    Loads sequences from a FASTA file into a dictionary: {contig_id: sequence}.
    """
    return {record.id: str(record.seq) for record in SeqIO.parse(fasta_path, "fasta")}

def load_aligned_contigs(snp_tsv):
    """
    Parses the snp_positions.tsv to extract pairs of (ref, alt) contigs.
    """
    df = pd.read_csv(snp_tsv, sep="\t")
    return list(set(zip(df["Target"], df["Query"])))  # Target = ref, Query = alt (check field order!)

def simulate_CO(ref_seq, alt_seq):
    """
    Simulates a single crossover between two parental sequences.
    """
    seq_len = len(ref_seq)
    crossover = np.random.binomial(1, 0.5)
    if crossover:
        crossover_pos = np.random.randint(0, seq_len)
        if np.random.binomial(1, 0.5):
            left, right = ref_seq[:crossover_pos], alt_seq[crossover_pos:]
            left_parent, right_parent = "ref", "alt"
        else:
            left, right = alt_seq[:crossover_pos], ref_seq[crossover_pos:]
            left_parent, right_parent = "alt", "ref"
        return left + right, (left_parent, right_parent, crossover_pos)
    else:
        if np.random.binomial(1, 0.5):
            return ref_seq, ("ref", "ref", None)
        else:
            return alt_seq, ("alt", "alt", None)

def simulate_f1_genome(ref_genome, alt_genome, contig_pairs, rep, output_dir):
    """
    Simulates a single F1 hybrid genome from matching contig pairs.
    """
    hybrid_record = []
    f1_table = []

    for ref_contig, alt_contig in contig_pairs:
        if ref_contig not in ref_genome or alt_contig not in alt_genome:
            print(f"Skipping: {ref_contig} or {alt_contig} not found in FASTA files")
            continue

        ref_seq = ref_genome[ref_contig]
        alt_seq = alt_genome[alt_contig]
        hybrid_seq, (left, right, pos) = simulate_CO(ref_seq, alt_seq)

        hybrid_record.append(
            SeqRecord(
                Seq(hybrid_seq),
                id=f"{ref_contig}_rep{rep}",
                description="F1_recombinant"
            )
        )
        f1_table.append({
            "rep": rep,
            "ref_chrom": ref_contig,
            "alt_chrom": alt_contig,
            "left": left,
            "right": right,
            "crossover_pos": pos if pos is not None else "None"
        })

    fasta_out = os.path.join(output_dir, f"f1_genome_rep{rep}.fna")
    with open(fasta_out, "w") as handle:
        SeqIO.write(hybrid_record, handle, "fasta-2line")

    return f1_table

def main(ref_fasta, alt_fasta, snp_tsv, output_dir, n_replicates):
    os.makedirs(output_dir, exist_ok=True)

    print("Loading parent FASTAs...")
    ref_genome = load_fasta(ref_fasta)
    alt_genome = load_fasta(alt_fasta)

    print("Parsing aligned contig pairs from SNP TSV...")
    contig_pairs = load_aligned_contigs(snp_tsv)

    all_f1_tables = []
    for rep in range(1, n_replicates + 1):
        print(f"Simulating F1 replicate {rep}...")
        f1_table = simulate_f1_genome(ref_genome, alt_genome, contig_pairs, rep, output_dir)
        all_f1_tables.extend(f1_table)

    df = pd.DataFrame(all_f1_tables)
    df.to_csv(os.path.join(output_dir, "f1_table.csv"), index=False)
    print(f"Done. F1 genomes written to: {output_dir}")

def parse_args():
    parser = argparse.ArgumentParser(description="Simulate recombined F1 hybrid genome(s) from parental FASTAs.")
    parser.add_argument("--ref-fasta", required=True, help="Path to parent 1 genome (.fna)")
    parser.add_argument("--alt-fasta", required=True, help="Path to parent 2 genome (.fna)")
    parser.add_argument("--snp-tsv", required=True, help="Path to SNP positions TSV from align_parents.py")
    parser.add_argument("--outdir", default="genomes/f1_simulations", help="Output directory for F1s")
    parser.add_argument("--reps", type=int, default=1, help="Number of replicate F1s to simulate")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    main(args.ref_fasta, args.alt_fasta, args.snp_tsv, args.outdir, args.reps)

#====
# EXAMPLE CLI
#====

#(base) riyarampalli@Riyas-MacBook-Pro-7 distortopia % python scripts/create_f1_hybrid.py \
  #--ref-dir input-data/Arabidopsis_thaliana/ncbi_dataset/data/GCA_000001735.2/GCA_000001735.2_TAIR10.1_genomic.fna \
  #--alt-dir input-data/Arabidopsis_lyrata/ncbi_dataset/data/GCA_000004255.1/GCA_000004255.1_v.1.0_genomic.fna \
  ##--snp-tsv genomes/par_alignments/snp_positions.tsv \
  ##--outdir genomes/f1_simulations \
  ##--reps 3

'''
Use to reduce scope for debugging. Allows user to avoid the error 
zsh: killed, which means process was forcefully terminated by the operating system 
— most often because it exceeded available memory (RAM). 
'''
#(base) riyarampalli@MacBookPro distortopia % head -n 10 genomes/par_alignments/snp_positions.tsv > genomes/par_alignments/snp_subset.tsv





