from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import argparse
import numpy as np
import pandas as pd

# Load the simulated F1 genome FASTA into a dictionary
def load_f1_fasta(fasta_path):
    return {record.id: str(record.seq) for record in SeqIO.parse(fasta_path, "fasta")}

# Simulate long reads and align to parents
def simulate_and_align(f1_genome, f1_table, rep, output_dir, read_len=15000, coverage=60):
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

    # Save simulated reads
    read_path = os.path.join(output_dir, f"f1_reads_rep{rep}.fasta")
    with open(read_path, "w") as handle:
        SeqIO.write(records, handle, "fasta-2line")

    # Align to reference and alt
    sam_ref = os.path.join(output_dir, f"rep{rep}_vs_ref.sam")
    sam_alt = os.path.join(output_dir, f"rep{rep}_vs_alt.sam")

    print(f"Mapping reads (rep {rep})...")
    os.system(f"minimap2 -ax map-pb {ref_fasta} {read_path} > {sam_ref}")
    os.system(f"minimap2 -ax map-pb {alt_fasta} {read_path} > {sam_alt}")

# Main entry point
def main(ref_fasta_arg, alt_fasta_arg, snp_tsv, output_dir, n_reps):
    global ref_fasta, alt_fasta
    ref_fasta = ref_fasta_arg
    alt_fasta = alt_fasta_arg
    os.makedirs(output_dir, exist_ok=True)
    f1_table = pd.read_csv(os.path.join(output_dir, "f1_table.csv"))

    for rep in range(1, n_reps + 1):
        print(f"Simulating long reads for rep {rep}...")
        fasta_path = os.path.join(output_dir, f"f1_genome_rep{rep}.fna")
        f1_genome = load_f1_fasta(fasta_path)
        simulate_and_align(f1_genome, f1_table, rep, output_dir)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref-fasta", required=True)
    parser.add_argument("--alt-fasta", required=True)
    parser.add_argument("--snp-tsv", required=True)
    parser.add_argument("--outdir", default="genomes/f1_simulations")
    parser.add_argument("--reps", type=int, default=1)
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

#Field            What it tells you                                   Source
##POS             Read alignment start on reference genome            SAM field 4
##crossover_pos   Breakpoint between parental segments in F1 genome   f1_table.csv

#(base) riyarampalli@MacBookPro distortopia % python scripts/long_reads.py \
  #--ref-fasta input-data/Arabidopsis_thaliana/ncbi_dataset/data/GCA_020911765.2/GCA_020911765.2_ASM2091176v2_genomic.fna \
  #--alt-fasta input-data/Arabidopsis_lyrata/ncbi_dataset/data/GCA_000524985.1/GCA_000524985.1_Alyr_1.0_genomic.fna \
  #--snp-tsv genomes/par_alignments/snp_positions.tsv \
  #--outdir genomes/f1_simulations \
  #--reps 1

'''
Print first few header lines of SAM file
'''
#head genomes/f1_simulations/rep1_vs_ref.sam

#@HD: Header line
##VN:1.6 = SAM format version 1.6
##SO:unsorted = reads are not sorted
##GO:query = grouped by query name

#@SQ: Sequence dictionary lines (one for each reference contig)
##SN:<name> = sequence name (e.g. CP087126.2)
##LN:<length> = contig length in bp (e.g. 32,638,291 bp)



