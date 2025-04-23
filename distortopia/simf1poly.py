import argparse
import random
from Bio import SeqIO
import pandas as pd

def load_fasta_dict(fasta_path):
    """
    Load sequences from FASTA into a dictionary: {contig: sequence string}
    """
    return {record.id.split()[0]: str(record.seq) for record in SeqIO.parse(fasta_path, "fasta")}

def parse_snp_table(tsv_path):
    """
    Read SNP positions table with columns: Query, Target, SNP_Positions
    Returns: list of tuples (query_contig, target_contig, [positions])
    """
    df = pd.read_csv(tsv_path, sep="\t")
    snp_entries = []
    for _, row in df.iterrows():
        if pd.isna(row["SNP_Positions"]) or str(row["SNP_Positions"]).strip() == "":
            continue
        positions = list(map(int, str(row["SNP_Positions"]).split(", ")))
        snp_entries.append((row["Query"], row["Target"], positions))
    return snp_entries

def simulate_f1(query_fasta, target_fasta, snp_entries, out_fasta, out_vcf):
    """
    Simulate F1 FASTA and VCF based on SNP positions across aligned contigs.
    """
    p1 = load_fasta_dict(query_fasta)   # Parent 1 = Query
    p2 = load_fasta_dict(target_fasta)  # Parent 2 = Target

    vcf_rows = []
    f1_contigs = {}

    for query_ctg, target_ctg, positions in snp_entries:
        if query_ctg not in p1 or target_ctg not in p2:
            print(f"[SKIP] Missing contig: {query_ctg} or {target_ctg}")
            continue

        seq1 = p1[query_ctg]
        seq2 = p2[target_ctg]
        min_len = min(len(seq1), len(seq2))

        f1_seq = list(seq1[:min_len])  # base: use query contig length

        for pos in positions:
            i = pos - 1  # convert to 0-based
            if i >= min_len:
                continue
            base1 = seq1[i]
            base2 = seq2[i]
            if base1 not in "ACGT" or base2 not in "ACGT":
                continue
            f1_base = random.choice([base1, base2])
            f1_seq[i] = f1_base

            if base1 != base2:
                vcf_rows.append([target_ctg, pos, ".", base1, base2, ".", ".", ".", "GT", f1_base])

        f1_contigs[query_ctg] = "".join(f1_seq)

    # Write pseudo-F1 FASTA
    with open(out_fasta, "w") as outfa:
        for contig, seq in f1_contigs.items():
            outfa.write(f">{contig}\n")
            for i in range(0, len(seq), 60):
                outfa.write(seq[i:i+60] + "\n")

    # Write VCF
    with open(out_vcf, "w") as outvcf:
        outvcf.write("##fileformat=VCFv4.2\n")
        outvcf.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Simulated haploid genotype\">\n")
        outvcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tF1\n")
        for row in vcf_rows:
            outvcf.write("\t".join(map(str, row)) + "\n")

    print(f"[✓] Simulated F1 FASTA written to: {out_fasta}")
    print(f"[✓] Simulated F1 VCF written to: {out_vcf}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulate a haploid F1 from aligned parental genomes using SNPs from a .tsv")
    parser.add_argument("--parent1", required=True, help="FASTA for species 1 (Query)")
    parser.add_argument("--parent2", required=True, help="FASTA for species 2 (Target)")
    parser.add_argument("--snps", required=True, help="TSV from summary_of_paf() with SNP positions")
    parser.add_argument("--out-fasta", default="f1_hybrid.fna", help="Output FASTA file for pseudo-F1")
    parser.add_argument("--out-vcf", default="f1_genotypes.vcf", help="Output VCF file of simulated F1 genotypes")
    args = parser.parse_args()

    snp_entries = parse_snp_table(args.snps)
    simulate_f1(args.parent1, args.parent2, snp_entries, args.out_fasta, args.out_vcf)

## ==========
# EXAMPLE CLI 
## ==========

"""
def run_simulation(args):
    parent1 = load_vcf(args.parent1) #this VCF generated from `simparents.py` script 
    parent2 = load_vcf(args.parent2) #this VCF generated from `simparents.py` script 

    f1_hybrid = {} #an empty dictionary? (key:value)

    for variant in shared_variants(parent1, parent2):
        allele1 = select_random_allele(parent1[variant])
        allele2 = select_random_allele(parent2[variant])
        f1_hybrid[variant] = f"{allele1}/{allele2}"

    write_vcf(f1_hybrid, output=args.output) #consider making this output a CSV?  

#propbabiltiy of croosover on each chormosome
"""