import argparse
import random
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

def generate_f1_hybrid(parent1_fasta, parent2_fasta, snp_zones, out_fasta, contig_map, recomb_per_contig=1):
    parent1 = parse_fasta_with_clean_ids(parent1_fasta)
    parent2 = parse_fasta_with_clean_ids(parent2_fasta)

    hybrid_records = []
    missing = 0
    for query_contig, target_contig in contig_map.items():
        if target_contig not in parent1:
            print(f"Skipping: {target_contig} not in parent1 FASTA")
            missing += 1
            continue
        if query_contig not in parent2:
            print(f"Skipping: {query_contig} not in parent2 FASTA")
            missing += 1
            continue

        seq1 = parent1[target_contig].seq
        seq2 = parent2[query_contig].seq

        zones = snp_zones.get(target_contig, [])
        breakpoints = []

        if zones:
            for _ in range(min(recomb_per_contig, len(zones))):
                bp = random.choice(zones)
                mid = (bp[0] + bp[1]) // 2
                breakpoints.append(mid)
        breakpoints = sorted(breakpoints)

        hybrid_seq = ""
        last = 0
        toggle = True
        for bp in breakpoints:
            hybrid_seq += seq1[last:bp] if toggle else seq2[last:bp]
            last = bp
            toggle = not toggle
        hybrid_seq += seq1[last:] if toggle else seq2[last:]

        hybrid_records.append(SeqIO.SeqRecord(Seq(hybrid_seq), id=f"{query_contig}_F1", description="Simulated F1 hybrid"))

    SeqIO.write(hybrid_records, out_fasta, "fasta")
    print(f"\nHybrid genome written to: {out_fasta}")
    print(f"{len(hybrid_records)} contigs written, {missing} skipped due to missing contigs")

def parse_args():
    parser = argparse.ArgumentParser(description="Simulate F1 hybrid using SNP-rich recombination zones")
    parser.add_argument("--ref-dir", required=True, help="path to parent 1 FASTA")
    parser.add_argument("--query-dir", required=True, help="path to parent 2 FASTA")
    parser.add_argument("--paf", required=True, help="path to .paf file output from minimap2 with cs:Z tags")
    parser.add_argument("--out", required=True, help="path to output hybrid FASTA")
    parser.add_argument("--recomb", type=int, default=1, help="number of recombination events per contig")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    print("Parsing .paf and extracting SNP zones...")
    snp_zones = extract_from_paf(args.paf)

    print("Mapping query contigs to target contigs...")
    contig_map = get_contig_mapping(args.paf)

    print("Building synthetic F1 hybrid...")
    generate_f1_hybrid(
        args.ref_dir,
        args.query_dir,
        snp_zones,
        args.out,
        contig_map,
        recomb_per_contig=args.recomb
    )

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