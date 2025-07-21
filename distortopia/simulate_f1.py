# distortopia/simulate_f1.py

import sys
from Bio import SeqIO

# IUPAC codes for heterozygous SNPs
IUPAC = {
    frozenset(['A', 'G']): 'R',
    frozenset(['C', 'T']): 'Y',
    frozenset(['G', 'C']): 'S',
    frozenset(['A', 'T']): 'W',
    frozenset(['G', 'T']): 'K',
    frozenset(['A', 'C']): 'M'
}

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

def write_fasta(reference, output_file):
    with open(output_file, "w") as out:
        for chrom in reference:
            seq = "".join(reference[chrom])
            out.write(f">{chrom}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")

def generate_f1_from_files(ref_fasta, vcf1, vcf2, output_path):
    ref_seq = load_reference(ref_fasta)
    variants1 = load_variants(vcf1)
    variants2 = load_variants(vcf2)
    f1_seq = apply_f1_variants(ref_seq, variants1, variants2)
    write_fasta(f1_seq, output_path)

# --- CLI mode ---
if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python simulate_f1.py ref.fna parent1.vcf parent2.vcf output.fna")
        sys.exit(1)
    generate_f1_from_files(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    print(f"F1_hybrid FASTA written to {sys.argv[4]}")

