import sys
from Bio import SeqIO

def load_reference(fasta_file):
    ref = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        ref[record.id] = list(record.seq)
    return ref

def apply_variants(reference, vcf_file):
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
    with open(output_file, "w") as out:
        for chrom in reference:
            seq = "".join(reference[chrom])
            out.write(f">{chrom}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i+60] + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python generate_f1_hybrid.py ref.fna variants.vcf output.fna")
        sys.exit(1)

    ref_fasta, variant_vcf, output_fasta = sys.argv[1], sys.argv[2], sys.argv[3]

    ref_seq = load_reference(ref_fasta)
    modified = apply_variants(ref_seq, variant_vcf)
    write_fasta(modified, output_fasta)
    print(f" F1 hybrid FASTA written to {output_fasta}")
