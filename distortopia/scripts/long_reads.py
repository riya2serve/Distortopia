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

def parse_sam_to_vcf(sam_path, ref_label, rep, f1_table):
    records = []
    for line in open(sam_path, "r"):
        if line.startswith("@"):
            continue
        fields = line.strip().split("\t")
        if len(fields) < 11:
            continue
        query_name = fields[0]
        chrom = fields[2]
        pos = fields[3]
        mapq = fields[4]

        # Extract contig name from read ID like "rep1_CP087129.2_read3"
        match = re.match(r"rep\d+_(.+?)_read\d+", query_name)
        if match:
            contig = match.group(1)
            row = f1_table[(f1_table["rep"] == rep) & (f1_table["ref_chrom"] == contig)]
            if not row.empty:
                left = row.iloc[0]["left"]
                right = row.iloc[0]["right"]
                crossover = row.iloc[0]["crossover_pos"]
                origin = f"{left}-{right}"
                xover_str = str(crossover) if crossover != "None" else "NA"
            else:
                origin = "unknown"
                xover_str = "NA"
        else:
            origin = "unknown"
            xover_str = "NA"

        records.append({
            "CHROM": chrom,
            "POS": pos,
            "ID": query_name,
            "REF": ".",  # placeholder
            "ALT": ".",  # placeholder
            "QUAL": mapq,
            "FILTER": "PASS",
            "INFO": f"source={ref_label};origin={origin};crossover_pos={xover_str}"
        })
    return records

def simulate_long_reads(f1_genome, f1_table, rep, output_dir, ref_fasta, alt_fasta, read_len=15000, coverage=60):
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

    # === Align to both parents ===
    sam_ref = os.path.join(output_dir, f"rep{rep}_vs_ref.sam")
    sam_alt = os.path.join(output_dir, f"rep{rep}_vs_alt.sam")

    print(f"Mapping long reads for rep{rep} to reference genomes...")
    os.system(f"minimap2 -ax map-pb {ref_fasta} {out_path} > {sam_ref}")
    os.system(f"minimap2 -ax map-pb {alt_fasta} {out_path} > {sam_alt}")

    # === Parse SAMs into unified VCF-like table ===
    vcf_like_records = parse_sam_to_vcf(sam_ref, "ref", rep, f1_table) + \
                       parse_sam_to_vcf(sam_alt, "alt", rep, f1_table)
    vcf_df = pd.DataFrame(vcf_like_records)
    vcf_out_path = os.path.join(output_dir, f"rep{rep}_long_read_alignment.vcf.tsv")
    vcf_df.to_csv(vcf_out_path, sep="\t", index=False)
    print(f"VCF-like alignment summary saved to: {vcf_out_path}")

def main(ref_fasta, alt_fasta, snp_tsv, output_dir, n_reps):
    os.makedirs(output_dir, exist_ok=True)
    f1_table = pd.read_csv(os.path.join(output_dir, "f1_table.csv"))

    for rep in range(1, n_reps + 1):
        print(f"Simulating long reads for rep {rep}...")
        fasta_path = os.path.join(output_dir, f"f1_genome_rep{rep}.fna")
        f1_genome = load_f1_fasta(fasta_path)
        simulate_long_reads(f1_genome, f1_table, rep, output_dir, ref_fasta, alt_fasta)

def parse_args():
    parser = argparse.ArgumentParser(description="Simulate and map long reads from F1 hybrid genomes.")
    parser.add_argument("--ref-fasta", required=True, help="Path to parent 1 FASTA")
    parser.add_argument("--alt-fasta", required=True, help="Path to parent 2 FASTA")
    parser.add_argument("--snp-tsv", required=True, help="Path to SNP table (unused but required for compatibility)")
    parser.add_argument("--outdir", default="genomes/f1_simulations", help="Directory to write outputs")
    parser.add_argument("--reps", type=int, default=1, help="Number of F1 replicates to simulate")
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

