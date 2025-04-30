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
    #loads query ('reference') and target FASTA sequences into dictionaries
    return {record.id: str(record.seq) for record in SeqIO.parse(fasta_path, "fasta")}

def simulate_CO (ref_seq, alt_seq):
    """
    Simualte a single crossover event between ref_seq and alt_seq.
    Return the recombinant sequence and a tuple (left_parent, right_parent, crossover_position).
    """
    seq_len = len(ref_seq) #getting length of the reference seq
    crossover = np.random.binomial(1, 0.5) #50% chance of a crossover event
    if crossover:
        crossover_pos = np.random.randint(0, seq_len) #position will occur somewhere btween start-end of contig
        if np.random.binomial(1, 0.5): #chance of crossover occuring
            left = ref_seq[:crossover_pos] #left of crossover will be ref_seq
            right = alt_seq[crossover_pos:] #right of crossover will be alt_seq
            left_parent, right_parent = "ref", "alt"
        else: #otherwise the below is actually the order
            left = alt_seq[:crossover_pos]
            right = ref_seq[crossover_pos:]
            left_parent, right_parent = "alt", "ref"
        return left + right, (left_parent, right_parent, crossover_pos)
    else:
        if np.random.binomial(1,0.5):
            return ref_seq, ("ref", "ref", None) #no crossover, ref_seq inherited
        else:
            return alt_seq, ("alt", "alt", None) #no crossover, alt_seq inherited

def simulate_f1_genome(ref_genome, alt_genome, rep, output_dir):
    """
    Generate a single F1 genome sequence.
    Return both its FASTA file (.fna) and a csv of inherited blocks.
    """

    hybrid_record = [] #initializing empty list of SeqRecord objects to write to FASTA file 
    f1_table = [] #initializing a list of dicts describing inheritance 

    #Iterate over contigs by index num
    for (ref_chrom, ref_seq), (alt_chrom, alt_seq) in zip(ref_genome.items(), alt_genome.items()):
        hybrid_seq, (left, right, pos) = simulate_CO(ref_seq, alt_seq)

        # Save recombinant sequence
        hybrid_record.append(
            SeqRecord(
                Seq(hybrid_seq.upper()), #joining generated FASTA sequences
                id=f"{ref_chrom}_rep{rep}",  # use reference contig name in ID
                description="F1_recombinant"
            )
        )
        f1_table.append({
            "rep": rep,
            "ref_chrom": ref_chrom,
            "alt_chrom": alt_chrom,
            "left": left,
            "right": right,
            "crossover_pos": pos if pos is not None else "None" #writing positions
        })

#Write out hybrid FASTA
    fasta_out = os.path.join(output_dir, f"f1_genome_rep{rep}.fna")
    with open(fasta_out, "w") as handle:
        SeqIO.write(hybrid_record, handle, "fasta-2line")

    return f1_table #printing table

def main(ref_dir, alt_dir, output_dir, n_replicates):
    os.makedirs(output_dir, exist_ok = True)

    print("Loading reference and alternate genomes . . .")
    ref_genome = load_fasta(ref_dir)
    alt_genome = load_fasta(alt_dir)

    all_f1_tables = [] #initializing list for replicate F1 genomes 
    for rep in range(1,n_replicates +1):
        print(f"Simulating F1 replicate {rep}...")
        f1_table = simulate_f1_genome(ref_genome, alt_genome, rep, output_dir)
        all_f1_tables.extend(f1_table)

    #Save full table
    df = pd.DataFrame(all_f1_tables)
    df.to_csv(os.path.join(output_dir, "f1_table.csv"), index = False)
    print("All F1 simulations complete. Outputs written to:", output_dir)

def parse_args():
    parser = argparse.ArgumentParser(description="Simulate recombined F1 hybrid FASTA(s) from two parents.")
    parser.add_argument("--ref-dir", required=True, help="Path to parent 1 genome (.fna)")
    parser.add_argument("--alt-dir", required=True, help="Path to parent 2 genome (.fna)")
    parser.add_argument("--outdir", default="genomes", help="Output directory for simulated genomes")
    parser.add_argument("--reps", type=int, default=1, help="Number of F1 replicates to generate")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    main(args.ref_dir, args.alt_dir, args.outdir, args.reps)

#====
# EXAMPLE CLI
#====

#(base) riyarampalli@Riyas-MacBook-Pro-7 distortopia % python scripts/create_f1_hybrid.py \
  #--ref-dir input-data/Arabidopsis_thaliana/ncbi_dataset/data/GCA_000001735.2/GCA_000001735.2_TAIR10.1_genomic.fna \
  #--alt-dir input-data/Arabidopsis_lyrata/ncbi_dataset/data/GCA_000004255.1/GCA_000004255.1_v.1.0_genomic.fna \
  #--outdir genomes/f1_simulations \
  #--reps 3




