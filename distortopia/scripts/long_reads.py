from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import os
import argparse
import numpy as np
import re

# Load the simulated F1 genome FASTA into a dictionary {contig_id: sequence}
def load_f1_fasta(fasta_path):
    return {record.id: str(record.seq) for record in SeqIO.parse(fasta_path, "fasta")}

# Simulate long reads from F1 hybrid genome
def simulate_long_reads(f1_genome, f1_table, rep, output_dir, read_len=15000, coverage=60):
    records = []
    subtable = f1_table[f1_table["rep"] == rep]  # Get rows corresponding to this replicate
    seen_contigs = set(subtable["ref_chrom"].unique())  # List of contigs to simulate reads from

    for contig_id in seen_contigs:
        hybrid_id = f"{contig_id}_rep{rep}"
        if hybrid_id not in f1_genome:
            continue
        seq = f1_genome[hybrid_id]
        seq_len = len(seq)
        num_reads = (seq_len * coverage) // read_len  # Number of reads based on coverage

        for i in range(num_reads):
            start = np.random.randint(0, seq_len - read_len + 1)  # Random start
            fragment = seq[start:start + read_len]  # Extract read
            read_id = f"rep{rep}_{contig_id}_read{i}"
            records.append(SeqRecord(Seq(fragment), id=read_id, description="simulated_long_read"))

    # Save simulated reads
    out_path = os.path.join(output_dir, f"f1_reads_rep{rep}.fasta")
    with open(out_path, "w") as out_handle:
        SeqIO.write(records, out_handle, "fasta-2line")

    # === Map reads to both parents ===
    sam_ref = os.path.join(output_dir, f"rep{rep}_vs_ref.sam")
    sam_alt = os.path.join(output_dir, f"rep{rep}_vs_alt.sam")

    print(f"Mapping reads (rep {rep}) to both reference genomes...")
    os.system(f"minimap2 -ax map-pb {ref_fasta} {out_path} > {sam_ref}")
    os.system(f"minimap2 -ax map-pb {alt_fasta} {out_path} > {sam_alt}")

    # === Parse SAM files to extract real SNPs using MD:Z tag ===
    def parse_sam_to_vcf(sam_path, ref_label):
        records = []
        with open(sam_path, "r") as sam:
            for line in sam:
                if line.startswith("@"):  # skip header lines
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 12:
                    continue

                query_name = fields[0]
                chrom = fields[2]
                pos = int(fields[3])
                mapq = fields[4]
                seq = fields[9]
                tags = fields[11:]

                # Extract MD:Z tag (which encodes mismatches)
                md_tag = next((tag for tag in tags if tag.startswith("MD:Z:")), None)
                if not md_tag:
                    continue
                md = md_tag.split(":")[-1]

                ref_pos = pos  # 1-based position from SAM
                i = 0  # index in read seq
                ref_base = None
                alt_base = None

                # Parse MD:Z tag to identify mismatches
                num = ""
                for char in md:
                    if char.isdigit():
                        num += char
                    else:
                        if num != "":
                            i += int(num)  # advance read index
                            ref_pos += int(num)  # advance reference pos
                            num = ""
                        if char.isalpha():
                            ref_base = char
                            if i < len(seq):
                                alt_base = seq[i]
                                # Save as VCF-like SNP
                                records.append({
                                    "CHROM": chrom,
                                    "POS": ref_pos,
                                    "ID": query_name,
                                    "REF": ref_base,
                                    "ALT": alt_base,
                                    "QUAL": mapq,
                                    "FILTER": "PASS",
                                    "INFO": f"source={ref_label}"
                                })
                            i += 1
                # If there's a trailing number, just ignore — it means perfect matches after last SNP

        return records

    # Collect records from both parental mappings
    vcf_like_records = parse_sam_to_vcf(sam_ref, "ref") + parse_sam_to_vcf(sam_alt, "alt")
    vcf_df = pd.DataFrame(vcf_like_records)

    # Save to VCF-like TSV
    vcf_out_path = os.path.join(output_dir, f"rep{rep}_long_read_alignment.vcf.tsv")
    vcf_df.to_csv(vcf_out_path, sep="\t", index=False)
    print(f"VCF-like alignment summary saved to: {vcf_out_path}")

# Main function
def main(ref_fasta_arg, alt_fasta_arg, snp_tsv, output_dir, n_reps):
    global ref_fasta, alt_fasta  # So they can be accessed in simulate_long_reads
    ref_fasta = ref_fasta_arg
    alt_fasta = alt_fasta_arg

    os.makedirs(output_dir, exist_ok=True)
    f1_table = pd.read_csv(os.path.join(output_dir, "f1_table.csv"))  # Contains contigs used in F1

    for rep in range(1, n_reps + 1):
        print(f"Simulating long reads for rep {rep}...")
        fasta_path = os.path.join(output_dir, f"f1_genome_rep{rep}.fna")
        f1_genome = load_f1_fasta(fasta_path)
        simulate_long_reads(f1_genome, f1_table, rep, output_dir)

# Parse command-line args
def parse_args():
    parser = argparse.ArgumentParser(description="Simulate long reads from F1 hybrid genomes and generate VCF-like SNP calls.")
    parser.add_argument("--ref-fasta", required=True, help="Path to parent 1 FASTA")
    parser.add_argument("--alt-fasta", required=True, help="Path to parent 2 FASTA")
    parser.add_argument("--snp-tsv", required=True, help="Path to SNP table (unused, kept for compatibility)")
    parser.add_argument("--outdir", default="genomes/f1_simulations", help="Output directory")
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

#Field            What it tells you                                   Source
##POS             Read alignment start on reference genome            SAM field 4
##crossover_pos   Breakpoint between parental segments in F1 genome   f1_table.csv

#(base) riyarampalli@MacBookPro distortopia % python scripts/long_reads.py \
  #--ref-fasta input-data/Arabidopsis_thaliana/ncbi_dataset/data/GCA_020911765.2/GCA_020911765.2_ASM2091176v2_genomic.fna \
  #--alt-fasta input-data/Arabidopsis_lyrata/ncbi_dataset/data/GCA_000524985.1/GCA_000524985.1_Alyr_1.0_genomic.fna \
  #--snp-tsv genomes/par_alignments/snp_positions.tsv \
  #--outdir genomes/f1_simulations \
  #--reps 1