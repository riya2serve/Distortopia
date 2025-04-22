import argparse
import random
from collections import defaultdict
from Bio import SeqIO

def extract_from_paf(paf_path, max_snp_gap=1000):
    """
    Parses .paf file for cs:Z tags to extract SNP positions. Then groups nearby SNPs into "recombination zones."
    Returns: dict of {contig: [(start, end), ...]}
    """
    snp_zones = defaultdict(list)

    with open(paf_path) as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 12:
                continue

            ref_contig = cols[5]
            ref_start = int(cols[7])

            cs_tag = ""
            for field in cols[12:]:
                if field.startswith("cs:Z:"):
                    cs_tag = field[5:]
                    break

            snps = [] #initialize list to collect parsed cs:Z tags 
            pos = ref_start #start at beginning of reference genome (in this case A. thaliana)
            i = 0 #start count at zero
            while i < len(cs_tag):
                if cs_tag[i] == ":": #exact base-pair match
                    i += 1 #add one to counter
                    num = ""
                    while i < len(cs_tag) and cs_tag[i].isdigit():
                        num += cs_tag[i]
                        i += 1 
                    pos += int(num)
                elif cs_tag[i] == "*": #single-nucleotide polymorphism
                    snps.append(pos)
                    pos += 1
                    i += 3
                elif cs_tag[i] == "-": #deletion
                    i += 1 #add 1 to counter
                    while i < len(cs_tag) and cs_tag[i].isalpha():
                        i += 1
                    pos += 1
                elif cs_tag[i] == "+": #insertion
                    i += 1 #add 1 to counter
                    while i < len(cs_tag) and cs_tag[i].isalpha():
                        i += 1
                else:
                    i += 1

            #group SNPs into "recombination zones"
            if snps:
                snps.sort()
                start = snps[0]
                end = start
                for s in snps[1:]:
                    if s - end <= max_snp_gap: #max_snp_gap is the ____
                        end = s
                    else:
                        snp_zones[ref_contig].append((start, end))
                        start = s
                        end = s
                snp_zones[ref_contig].append((start, end))

    return snp_zones

def generate_f1_hybrid(parent1_fasta, parent2_fasta, snp_zones, out_fasta, recomb_per_contig=1):
    """
    Simulates recombination between two parental genomes based on SNP-rich zones.
    """
    parent1 = SeqIO.to_dict(SeqIO.parse(parent1_fasta, "fasta"))
    parent2 = SeqIO.to_dict(SeqIO.parse(parent2_fasta, "fasta"))

    hybrid_records = []

    for contig in parent1:
        if contig not in parent2:
            continue

        seq1 = parent1[contig].seq
        seq2 = parent2[contig].seq

        zones = snp_zones.get(contig, [])
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

        hybrid_records.append(SeqIO.SeqRecord(hybrid_seq, id=f"{contig}_F1", description="Simulated F1 hybrid"))

    SeqIO.write(hybrid_records, out_fasta, "fasta")
    print(f"Hybrid genome written to: {out_fasta}")

def parse_args():
    parser = argparse.ArgumentParser(description="Simulate F1 hybrid using SNP-rich recombination zones")
    parser.add_argument("--parent1", required=True, help="path to parent 1 FASTA")
    parser.add_argument("--parent2", required=True, help="path to parent 2 FASTA")
    parser.add_argument("--paf", required=True, help="path to .paf file output from minimap2 with cs:Z tags")
    parser.add_argument("--out", required=True, help="path to output hybrid FASTA")
    parser.add_argument("--recomb", type=int, default=1, help="Number of recombination events per contig")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    print("Parsing .paf and extracting SNP zones...")
    snp_zones = extract_from_paf(args.paf)

    print("Building synthetic F1 hybrid...")
    generate_f1_hybrid(args.parent1, args.parent2, snp_zones, args.out, recomb_per_contig=args.recomb)

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