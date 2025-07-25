# distortopia/simulate_f1.py

import sys
from Bio import SeqIO
from itertools import chain
import numpy as np
import bisect


def load_reference(fasta_file):
    ref = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        ref[record.id] = list(record.seq.upper())
    return ref

def load_variants(vcf_file):
    variants = {}
    with open(vcf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 5:
                continue
            chrom, pos, _, ref, alt = fields[:5]
            pos = int(pos) - 1  # VCF is 1-based
            if len(ref) == 1 and len(alt) == 1:
                variants.setdefault(chrom, {})[pos] = alt.upper()
    return variants

def apply_f1_variants(reference, var1, var2):
    ref1 = reference
    ref2 = reference.copy()

    # ...
    for chrom in reference:
        if chrom not in var1 and chrom not in var2:
            continue
        for pos in range(len(reference[chrom])):
            ref_base = reference[chrom][pos]
            base1 = var1.get(chrom, {}).get(pos, ref_base)
            base2 = var2.get(chrom, {}).get(pos, ref_base)

            if base1 == base2:
                reference[chrom][pos] = base1 #homozygous; use this base
            elif base1 == ref_base:
                reference[chrom][pos] = base2 #one parent differs; use ALT base
            elif base2 == ref_base:
                reference[chrom][pos] = base1
            else:
                code = IUPAC.get(frozenset([base1, base2]), 'N') #heterozygous; so use IUPAC
                reference[chrom][pos] = code
    return reference


def crossover_random(reference, var1, var2):
    ref = reference.copy()
    for chrom in ref:
        chrom_len = len(chrom)        
        # random sample whether a crossover occurs on this chrom
        if not np.random.binomial(1, 0.5):
            # randomly no crossover pos as 0 or end of chrom pos
            if np.random.binomial(1, 0.5):
                crossover_pos = 0
            else:
                crossover_pos = chrom_len
        else:
            # random uniform sample position of crossover
            crossover_pos = np.random.uniform(0, chrom_len)
            
        # apply variants from var1 until pos, then apply variants from var2
        for pos in range(len(ref[chrom])):
            ref_base = ref[chrom][pos]
            if pos < crossover_pos:
                base = var1.get(chrom, {}).get(pos, ref_base)
            else:
                base = var2.get(chrom, {}).get(pos, ref_base)
            ref[chrom][pos] = base
    return ref


def crossover_random_and_sim_reads(reference, var1, var2):
    ref = reference.copy()
    for chrom in ref:
        chrom_len = len(chrom)        
        # random sample whether a crossover occurs on this chrom
        if not np.random.binomial(1, 0.5):
            # randomly no crossover pos as 0 or end of chrom pos
            if np.random.binomial(1, 0.5):
                crossover_pos = 0
            else:
                crossover_pos = chrom_len
        else:
            # random uniform sample position of crossover
            crossover_pos = np.random.uniform(0, chrom_len)
            
        # apply variants from var1 until pos, then apply variants from var2
        for pos in range(len(ref[chrom])):
            ref_base = ref[chrom][pos]
            if pos < crossover_pos:
                base = var1.get(chrom, {}).get(pos, ref_base)
            else:
                base = var2.get(chrom, {}).get(pos, ref_base)
            ref[chrom][pos] = base

    # sample/sim reads from this reference
    reads = []
    chrom_lens = [len(i) for i in ref]

    # iterate over copies of the genome
    # get list of positions to sim reads from
    # assuming 1000 reps (genomes)
    # and readlen=10000
    # we need 3000 reads per rep to get 60X coverage.
    for i in range(30):
        chrom_idx, pos_in_chrom = sample_chrom_position(chrom_lens)
        try:
            chrom = ref[chrom_idx].copy()
            read = chrom[pos_in_chrom:pos_in_chrom + 10000]
            read.id = f"read{i}"
            reads.append(read)
        except IndexError:
            pass
    return reads


def sample_chrom_position(chrom_lens):
    chrom_starts = np.cumsum([0] + chrom_lens[:-1]).tolist()
    total_len = sum(chrom_lens)
    print(f"Total length:{chrom_lens}")
    rand_pos = np.random.randint(0, total_len)
    chrom_index = bisect.bisect_right(chrom_starts, rand_pos) - 1
    pos_in_chrom = rand_pos - chrom_starts[chrom_index]
    return chrom_index, pos_in_chrom


def write_fasta(copies, output_file):
    with open(output_file, "w") as out:
        for rep, reference in copies.items():
            for chrom in reference:
                seq = reference[chrom]
                seqstr = "".join(seq[chrom])
                out.write(f">{chrom}_rep{rep}_{seq.id}\n")
                for i in range(0, len(seqstr), 60):
                    out.write(seqstr[i:i+60] + "\n")

#def generate_f1_from_files(ref_fasta, vcf1, vcf2, output_path):
    #ref_seq = load_reference(ref_fasta)
    #variants1 = load_variants(vcf1)
    #variants2 = load_variants(vcf2)
    #f1_seq = apply_f1_variants(ref_seq, variants1, variants2)
    #write_fasta(f1_seq, output_path)

def generate_f1_from_files(ref_fasta, vcf1, vcf2, output_path, num_reps):
    ref_seq = load_reference(ref_fasta)
    variants1 = load_variants(vcf1)
    variants2 = load_variants(vcf2)
    copies = {}
    for i in range(num_reps): #generating multiple F1 fastas
        f1_seq = crossover_random_and_sim_reads(ref_seq, variants1, variants2)
        copies[i] = f1_seq
    write_fasta(copies, output_path)

# --- CLI mode ---
if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python simulate_f1.py ref.fna parent1.vcf parent2.vcf output.fna")
        sys.exit(1)
    generate_f1_from_files(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    print(f"F1_hybrid FASTA written to {sys.argv[4]}")

