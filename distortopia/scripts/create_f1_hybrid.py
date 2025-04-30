from Bio import SeqIO #BioPython module for reading and writing sequence files 
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord
import pandas as pd #for working with dataframes and generated HTML style summaries
import os #interacts with operating system 
import argparse #allows for command-line arguments/flags from/for users
import numpy as np #for assigning random binomial

"""
This script loads A. thaliana and A. lyrata .fna files and then parses the SNP alignment summary.
It allows users to extract both reference and target bases at SNP positions. 
With this, users can simulate a haploid F1 hybrid, with target SNPs incorporated into the reference
genome. 
"""
def load_fasta(fasta_path):
    # Load sequences into dictionary: {contig_id: sequence string}
    return {record.id: str(record.seq) for record in SeqIO.parse(fasta_path, "fasta")}

def load_aligned_contigs(snp_tsv):
    """
    Reads aligned contig pairs from the summary TSV file.
    Returns a list of (ref_contig, alt_contig) tuples.
    """
    df = pd.read_csv(snp_tsv, sep="\t")
    return list(zip(df["Query"], df["Target"]))

def simulate_CO(ref_seq, alt_seq): #for crossovers
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
    hybrid_record = []
    f1_table = []

    for ref_chrom, alt_chrom in contig_pairs:
        if ref_chrom not in ref_genome or alt_chrom not in alt_genome:
            continue
        ref_seq = ref_genome[ref_chrom]
        alt_seq = alt_genome[alt_chrom]

        hybrid_seq, (left, right, pos) = simulate_CO(ref_seq, alt_seq)

        hybrid_record.append(
            SeqRecord(
                Seq(hybrid_seq.upper()),
                id=f"{ref_chrom}_rep{rep}",
                description="F1_recombinant"
            )
        )
        f1_table.append({
            "rep": rep,
            "ref_chrom": ref_chrom,
            "alt_chrom": alt_chrom,
            "left": left,
            "right": right,
            "crossover_pos": pos if pos is not None else "None"
        })

    fasta_out = os.path.join(output_dir, f"f1_genome_rep{rep}.fna")
    with open(fasta_out, "w") as handle:
        SeqIO.write(hybrid_record, handle, "fasta-2line")

    return f1_table

def main(ref_dir, alt_dir, snp_tsv, output_dir, n_replicates):
    os.makedirs(output_dir, exist_ok=True)

    print("Loading reference and alternate genomes...")
    ref_genome = load_fasta(ref_dir)
    alt_genome = load_fasta(alt_dir)
    contig_pairs = load_aligned_contigs(snp_tsv)

    all_f1_tables = []
    for rep in range(1, n_replicates + 1):
        print(f"Simulating F1 replicate {rep}...")
        f1_table = simulate_f1_genome(ref_genome, alt_genome, contig_pairs, rep, output_dir)
        all_f1_tables.extend(f1_table)

    df = pd.DataFrame(all_f1_tables)
    df.to_csv(os.path.join(output_dir, "f1_table.csv"), index=False)
    print("All F1 simulations complete. Outputs written to:", output_dir)

def parse_args():
    parser = argparse.ArgumentParser(description="Simulate F1 hybrids using aligned contigs from SNP TSV.")
    parser.add_argument("--ref-dir", required=True, help="Path to parent 1 genome (.fna)")
    parser.add_argument("--alt-dir", required=True, help="Path to parent 2 genome (.fna)")
    parser.add_argument("--snp-tsv", required=True, help="Path to snp_positions.tsv from align_parents.py")
    parser.add_argument("--outdir", default="genomes", help="Output directory")
    parser.add_argument("--reps", type=int, default=1, help="Number of F1 replicates to generate")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    main(args.ref_dir, args.alt_dir, args.snp_tsv, args.outdir, args.reps)
#====
# EXAMPLE CLI
#====

#(base) riyarampalli@Riyas-MacBook-Pro-7 distortopia % python scripts/create_f1_hybrid.py \
  #--ref-dir input-data/Arabidopsis_thaliana/ncbi_dataset/data/GCA_000001735.2/GCA_000001735.2_TAIR10.1_genomic.fna \
  #--alt-dir input-data/Arabidopsis_lyrata/ncbi_dataset/data/GCA_000004255.1/GCA_000004255.1_v.1.0_genomic.fna \
  ##--snp-tsv genomes/par_alignments/snp_positions.tsv \
  ##--outdir genomes/f1_simulations \
  ##--reps 3


