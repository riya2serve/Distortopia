import sys
from Bio import SeqIO

def load_reference(fasta_file):
    """Loads a reference genome into a dictionary."""
    ref = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        ref[record.id] = list(record.seq)
    return ref

def apply_variants(reference, vcf_file):
    "Applies SNP variants identified in VCF file to reference genome."""
    with open(vcf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            chrom, pos, _, ref, alt, *_ = line.strip().split("\t")
            pos = int(pos) - 1  # VCF is 1-based
            if chrom in reference and reference[chrom][pos] == ref:
                reference[chrom][pos] = alt
    return reference

def write_fasta(reference, output_file):
    """Writes the modified reference genome to a new FASTA file."""
    with open(output_file, "w") as out:
        for chrom in reference:
            seq = "".join(reference[chrom])
            out.write(f">{chrom}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")

def generate_f1_hybrid(ref_fasta, vcf_file, output_fasta):
    """Callable function to simulate the hybrids."""
    ref_seq = load_reference(ref_fasta)
    modified = apply_variants(ref_seq, variant_vcf)
    write_fasta(modified, output_fasta)
    print(f"F1 hybrid FASTA written to {output_fasta}")
