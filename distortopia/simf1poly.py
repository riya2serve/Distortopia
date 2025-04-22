import argparse #allows for command-line arguments/flags from/for users
import random #to simulate recombination events 
from collections import defaultdict #to store SNP count per contig
from Bio import SeqIO #BioPython module for reading and writing sequence files 

def extract_from_paf(paf_path, max_snp_gap=1000):
    """
    Parses .paf file for cs:Z tags to extract SNP positions. Then groups nearby SNPs into "recombination zones."
    Returns: dict of {contig: [(start, end), ...]}
    """
    snp_zones = defaultdict(list) #initializing dict to store SNP regions by contig 

    with open(paf_path) as f: #opens .paf file (should already be generated)
        for line in f:
            cols = line.strip().split("\t") #split .paf by tabs
            if len(cols) < 12: 
                continue #skip lines without expected fields

            ref_contig = cols[5] #contig name from reference genome 
            ref_start = int(cols[7]) #start position on reference

            #extracts cs:Z tag -- encodes alignment differences (SNPs, indels)
            cs_tag = ""
            for field in cols[12:]:
                if field.startswith("cs:Z:"): 
                    cs_tag = field[5:] #remove cs:Z line prefix
                    break

            snps = [] #initialize list to collect SNP positions 
            pos = ref_start #track current position on reference genome (in this case A. thaliana)
            i = 0 #position in cs tag string
            while i < len(cs_tag):
                if cs_tag[i] == ":": #exact base-pair match
                    i += 1 #add one to counter
                    num = ""
                    while i < len(cs_tag) and cs_tag[i].isdigit(): #read full number of matched based 
                        num += cs_tag[i]
                        i += 1 
                    pos += int(num) #advance position on reference
                elif cs_tag[i] == "*": #single-nucleotide polymorphism
                    snps.append(pos) #record SNP position 
                    pos += 1 #advance position on reference
                    i += 3 #skip ref + alt bases 
                elif cs_tag[i] == "-": #deletion
                    i += 1 #add 1 to counter
                    while i < len(cs_tag) and cs_tag[i].isalpha(): #skip deleted bases 
                        i += 1
                    pos += 1 #advance position on reference 
                elif cs_tag[i] == "+": #insertion
                    i += 1 #add 1 to counter
                    while i < len(cs_tag) and cs_tag[i].isalpha(): #skip inserted bases
                        i += 1
                else:
                    i += 1 #catch all

            #group nearby SNPs into "recombination zones"
            if snps:
                snps.sort()
                start = snps[0]
                end = start
                for s in snps[1:]:
                    if s - end <= max_snp_gap: #if the SNP is close enough, then extend "recombination zone"
                        end = s
                    else:
                        snp_zones[ref_contig].append((start, end)) #save "recombination zone"
                        start = s
                        end = s
                snp_zones[ref_contig].append((start, end)) #add zones

    return snp_zones #return dictionary of "recombination zones" per contig
    """results.append({
                "Query": ...
    """

def generate_f1_hybrid(parent1_fasta, parent2_fasta, snp_zones, out_fasta, recomb_per_contig=1):
    """
    Simulates recombination between two parental genomes based on SNP-rich zones.
    Writes a new hybrid FASTA file.
    """
    parent1 = SeqIO.to_dict(SeqIO.parse(parent1_fasta, "fasta"))
    parent2 = SeqIO.to_dict(SeqIO.parse(parent2_fasta, "fasta"))

    hybrid_records = [] #initializing list to store hybrid contigs 

    for contig in parent1:
        if contig not in parent2:
            continue #skip contigs that are not shared by both parents 

        seq1 = parent1[contig].seq #parent 1 sequence
        seq2 = parent2[contig].seq #parent 2 sequence 

        zones = snp_zones.get(contig, []) #fetch SNPs zones for each contig
        breakpoints = []

        #choose recombination breakpoints randomly based on SNP-rich zones
        if zones:
            for _ in range(min(recomb_per_contig, len(zones))):
                bp = random.choice(zones)
                mid = (bp[0] + bp[1]) // 2 #use midpoint of "recombination zone" as breakpoint
                breakpoints.append(mid)
        breakpoints = sorted(breakpoints) #ordered list of breakpoints 

        #build hybrid sequence by recombining segments between the two parents
        hybrid_seq = ""
        last = 0
        toggle = True #determines which parent 'chunk' will be contributed
        for bp in breakpoints:
            hybrid_seq += seq1[last:bp] if toggle else seq2[last:bp]
            last = bp
            toggle = not toggle #alternated parent 'chunk'
        hybrid_seq += seq1[last:] if toggle else seq2[last:] #add final sequence 

        hybrid_records.append(SeqIO.SeqRecord(hybrid_seq, id=f"{contig}_F1", description="Simulated F1 hybrid"))

    #write all hybrid contigs to an output FASTA
    SeqIO.write(hybrid_records, out_fasta, "fasta")
    print(f"Hybrid genome written to: {out_fasta}")

def parse_args():
    """
    Command-line argument interface.
    Users can provide paths, number of threads, and minimap2 settings.
    """
    parser = argparse.ArgumentParser(description="Simulate F1 hybrid using SNP-rich recombination zones")
    parser.add_argument("--ref-dir", required=True, help="path to parent 1 FASTA")
    parser.add_argument("--query-dir", required=True, help="path to parent 2 FASTA")
    parser.add_argument("--paf", required=True, help="path to .paf file output from minimap2 with cs:Z tags")
    parser.add_argument("--out", required=True, help="path to output hybrid FASTA")
    parser.add_argument("--recomb", type=int, default=1, help="number of recombination events per contig")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    #Step 1. extract SNP regions 
    print("Parsing .paf and extracting SNP zones...")
    snp_zones = extract_from_paf(args.paf)

    print("Building synthetic F1 hybrid...")
    generate_f1_hybrid(args.ref_dir, args.query_dir, snp_zones, args.out, recomb_per_contig=args.recomb)

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