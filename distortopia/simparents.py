import os
import subprocess
import argparse

"""
This script will take two input FASTA files (--ref and --query).
It will align the query to the reference using 'minimap2.'
It will then call variants using paftools.js
Variant calls will be ouput onto a real VCF, which can be used later
to simulate F1 hybrid genotypes or detect segregation distortion."
"""

def align_call_variants(ref_fasta, query_fasta, output_vcf):
    """
    Align query genome to reference and call variants.
    """
    paf_file = output_vcf.replace(".vcf", ".paf")

    print(f"Aligning {query_fasta} to {ref_fasta}")
    with open(paf_file, "w") as paf_out:
        subprocess.run([
            "minimap2", "-x", "asm5", ref_fasta, query_fasta
        ], stdout=paf_out, check=True)

    print(f"Calling SNPs using paftools.js")
    with open(output_vcf, "w") as vcf_out:
        subprocess.run([
            "paftools.js", "call", paf_file
        ], stdout=vcf_out, check=True)

    print(f"SNPs written to: {output_vcf}")


def parse_args():
    parser = argparse.ArgumentParser(description="Call real SNPs between two genomes using minimap2 + paftools.js")
    parser.add_argument("--ref", required=True, help="Reference genome FASTA file")
    parser.add_argument("--query", required=True, help="Query genome FASTA file")
    parser.add_argument("--vcf", required=True, help="Path to output VCF file")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    if not os.path.isfile(args.ref):
        print(f"Reference genome not found: {args.ref}")
        exit(1)
    if not os.path.isfile(args.query):
        print(f"Query genome not found: {args.query}")
        exit(1)

    os.makedirs(os.path.dirname(args.vcf), exist_ok=True)
    align_call_variants(args.ref, args.query, args.vcf)

