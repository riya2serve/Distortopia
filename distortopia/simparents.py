from Bio import Entrez, SeqIO #BioPython for NCBI access and FASTA parsing
import os 
import random #for random number gneration (will use for SNP positions)
import pandas as pd #for dataframe handling (not necessary)
import argparse

# ================
# CONFIGURATION
# ================

"""
This script will allow users to simulate two parental genomes from reference FASTA files.
It will introduce random SNPs to generate Parent 1 and Parent 2.

With these, the user can then simulate an F1 hybrid VCF. If SNPs differ, the F1 gets heterozygous (0/1);
if SNPs are identical the F1 gets homozygous (1/1 OR 0/0).
"""

genome_fasta = "spinach_genome/spinach_wg.fasta"
#consider providing 2 reference fastas for cases where samples are v divergent?
#divergent cases might have different reference files

output_folder = "spinach_genome" #edit this
os.makedirs(output_folder, exist_ok = True) #checks that output folder exists

# ================
# FUNCTION(s)
# ================

def simulate_snps(fasta_file, snp_count, output_vcf):
    """
    Simulate random SNPs in a genome FASTA and save as a VCF.
    """
    sequences = list(SeqIO.parse(fasta_file, "fasta")) #can contain contigs or chromosomes
    tot_nome_length = sum(len(seq) for seq in sequences) #add up the contig or chromosome length so you can know full range for 'placing' SNPs

    print(f"[INFO] Loaded genome: {tot_nome_length:,} bp across {len(sequences)} chromosomes")

    # Generate random SNP positions
    snp_positions = random.sample(range(1, tot_nome_length), snp_count) #randomly sample snps across entire ___ length
    snp_positions.sort() #sort snp positions 

    # Store SNPs per chromosome
    snps_on_chrom = {} #this is a dictionary for storing the SNPs on each chromsome
    genome_cursor = 0 #position on genome; initializing at zero

    for seq in sequences: #this for-loop iterates through each chromosome
        chrom = seq.id #name of the chromosome
        chrom_length = len(seq)
        chrom_snps = [] #an empty list to store the SNPs found on a chromsome

        #will keep looping as long as these conditions are true:
        while snp_positions and genome_cursor < snp_positions[0] <= genome_cursor + chrom_length:
            pos_in_chrom = snp_positions.pop(0) - genome_cursor #converting genome wide SNPs into chrom-relative position
            ref_base = seq.seq[pos_in_chrom - 1].upper() #gets the reference nucleotide and converts it to uppercase

            if ref_base not in "ACGT":
                continue  # Skip ambiguous/unknown bases

            # Generate a random alternative allele
            alt_base = random.choice([b for b in "ACGT" if b != ref_base])
            chrom_snps.append((pos_in_chrom, ref_base, alt_base)) #add a tuple describing the SNP

        genome_cursor += chrom_length

        if chrom_snps:
            snps_on_chrom[chrom] = chrom_snps #storing SNPS for each chromosome
            ###{
               ##"NC_079487.1": [(pos1, ref, alt), (pos2, ref, alt)],
               ##"NC_079488.1": [(pos3, ref, alt), ...]
            ###}

    # Write VCF file
    with open(output_vcf, "w") as vcf:
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write(f"##source=Simulated_Parent\n")
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPARENT1\tPARENT2\n")

    #want to make parents polymorphic 
        genotypes = [("0/0", "1/1"), ("1/1", "0/0"), ("0/1", "0/1")]  #genotype options
        #both homozygous OR both heterozygous
        
        for chrom, snps in snps_on_chrom.items():
            for pos, ref, alt in snps:
                gt_parent1, gt_parent2 = random.choice(genotypes)

                vcf.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt_parent1}\t{gt_parent2}\n") #heterozygous parents

    print(f"[DONE] Simulated SNPs written to {output_vcf}")

# ================
# ARGPARSE
# ================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Simulate parental variants by generating SNPs from a ref FASTA genome.")

    # Required args
    parser.add_argument('--input-fasta', type=str, required=True,
                        help='Path to the reference FASTA genome file')
    parser.add_argument('--parent1', type=str, required=True,
                        help='Output VCF path for Parent 1')
    parser.add_argument('--parent2', type=str, required=True,
                        help='Output VCF path for Parent 2')

    # Optional args
    parser.add_argument('--snp-count', type=int, default=1000,
                        help='Number of SNPs per parent (default: 1000)')

    return parser.parse_args()

# ================
# FUNC. EXECUTION
# ================
if __name__ == "__main__":
    num_snps = 1000
    args = parse_args()

    # Check input FASTA file exists
    if not os.path.isfile(args.input_fasta):
        print(f"[ERROR] Input FASTA file not found: {args.input_fasta}")
        exit(1)
    # Create output folder if necessary (based on Parent1 VCF path)
    output_folder = os.path.dirname(args.parent1)
    os.makedirs(output_folder, exist_ok=True)

    print("[INFO] Starting parental genome simulation...")

    # Simulate SNPs for Parent1
    simulate_snps(
        fasta_file=genome_fasta,
        snp_count=num_snps,
        output_vcf=os.path.join(output_folder, "parent1.vcf")
    )

    # Simulate SNPs for Parent2 (different SNP set)
    simulate_snps(
        fasta_file=genome_fasta,
        snp_count=num_snps,
        output_vcf=os.path.join(output_folder, "parent2.vcf")
    )

    print("[COMPLETE] Parental VCF simulation finished!")

"""
===========================================================
CHECK OUTPUT IN TERMINAL (POST-RUN)
===========================================================

After running this script and generating the parental VCF files, you can 
inspect and verify them directly from your terminal. Note: GitHubs file size limit is 100.00Mb. 

If your FASTA file size >100mb you will need either Git LFS, a git extension
for working with large files, or .gitignore, a text file that tells Git to ignore 
certain files or directories. Both are useful for preventing temporary files, large files, or 
sensitive data from being added to your repo. They keep your Git history lightweight!

1. View the first few lines of Parent1 VCF:
-------------------------------------------------
$ head spinach_genome/parent1.vcf #prints the top 10 lines 
GT = genotypes

2. Scroll through the file interactively:
-------------------------------------------------
$ less spinach_genome/parent1.vcf #prints >10 lines
(press 'q' to quit out of `less`)

3. Count the number of simulated SNPs (excluding headers):
-------------------------------------------------
$ grep -v "^#" spinach_genome/parent1.vcf | wc -l

4. Search for variants on a specific chromosome (example):
-------------------------------------------------
$ grep "NC_079487.1" spinach_genome/parent1.vcf

===========================================================
"""
