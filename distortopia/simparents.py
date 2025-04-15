import os
import glob
from Bio import SeqIO
import argparse
import pandas as pd #to genrate HTML with color-codes SNPs

"""
This script allows a user to select one genome FASTA from each of two species,
then compares them base-by-base for all matching contigs (by order/length).

It outputs a real, biologically meaningful VCF of SNPs between the two genome assemblies.
"""
def choose_fasta(species_name, species_dir):
    fasta_paths = glob.glob(f"{species_dir}/**/*genomic.fna", recursive=True)
    if not fasta_paths:
        raise FileNotFoundError(f"No FASTA files found for {species_name} in {species_dir}")

    print(f"\n[SELECT FASTA FOR {species_name}]")
    for i, path in enumerate(fasta_paths):
        print(f"{i+1}. {path}")

    while True:
        try:
            index = int(input(f"Enter number [1-{len(fasta_paths)}]: ")) - 1
            if 0 <= index < len(fasta_paths):
                return fasta_paths[index]
            else:
                print("Invalid selection. Try again.")
        except ValueError:
            print("Please enter a valid number.")

def compare_contigs(ref_fasta, query_fasta, output_vcf, mode="length", top_n=5):
    """
    Compare top N contigs by order or length between two genome assemblies.
    """
    ref_seqs = list(SeqIO.parse(ref_fasta, "fasta"))
    query_seqs = list(SeqIO.parse(query_fasta, "fasta"))

    if mode == "length":
        ref_seqs.sort(key=lambda s: -len(s))
        query_seqs.sort(key=lambda s: -len(s))

    pair_count = min(top_n, len(ref_seqs), len(query_seqs))
    print(f"\nComparing top {pair_count} contigs by {mode}...\n")

    with open(output_vcf, "w") as vcf:
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write("##source=BiopythonSNPComparer\n")
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for i in range(pair_count):
            ref = ref_seqs[i]
            query = query_seqs[i]
            chrom = f"{ref.id}_vs_{query.id}"

            ref_seq = str(ref.seq).upper()
            query_seq = str(query.seq).upper()
            min_len = min(len(ref_seq), len(query_seq))
            snp_count = 0

            for j in range(min_len):
                r = ref_seq[j]
                q = query_seq[j]
                if r != q and r in "ACGT" and q in "ACGT":
                    vcf.write(f"{chrom}\t{j+1}\t.\t{r}\t{q}\t.\tPASS\t.\n")
                    snp_count += 1

            print(f"Compared {chrom} → {snp_count:,} SNPs")

    print(f"\nSNP comparison complete → {output_vcf}")

def gen_HTML(vcf_path, html_out = "snp_summary.hmtl"):
    """
    Loads a user-generated VCF file, summarizes SNP per contig, exports a styled HTML table.
    """ 
    with open(vcf_path) as f:
        lines = [line for line in f if not line.startswith("#")]

    df = pd.read_csv(
        pd.compat.StringIO("".join(lines)),
        sep="\t", #separating columns 
        names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"] #column names
    )

    # Filter biallelic SNPs
    valid = {"A", "C", "G", "T"} #looking for base pairs
    df = df[(df["REF"].isin(valid)) & (df["ALT"].isin(valid)) & (df["REF"] != df["ALT"])]

    # Count SNPs per contig
    summary = df["CHROM"].value_counts().reset_index()
    summary.columns = ["CHROM", "SNP_COUNT"]

    # Color styling
    def color(val):
        if val > 50000:
            return "background-color: #d73027; color: white"
        elif val > 10000:
            return "background-color: #fc8d59"
        elif val > 1000:
            return "background-color: #fee08b"
        else:
            return "background-color: #d9ef8b"

    styled = summary.style.applymap(color, subset=["SNP_COUNT"])
    styled.to_html(html_out, index=False)

    print(f"SNP summary table written to: {html_out}")

def parse_args():
    parser = argparse.ArgumentParser(description="Compare top contigs from two genomes by length or order")
    parser.add_argument("--ref-dir", required=True, help="Directory containing reference FASTAs")
    parser.add_argument("--query-dir", required=True, help="Directory containing query FASTAs")
    parser.add_argument("--out", required=True, help="Output VCF file path")
    parser.add_argument("--mode", choices=["length", "order"], default="length",
                        help="How to match contigs: by length (default) or order")
    parser.add_argument("--top-n", type=int, default=5, help="Number of contigs to compare (default: 5)")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    ref_fasta = choose_fasta("Species 1", args.ref_dir)
    query_fasta = choose_fasta("Species 2", args.query_dir)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    compare_contigs(ref_fasta, query_fasta, args.out, mode=args.mode, top_n=args.top_n)
    gen_HTML(args.out, html_out=args.out.replace(".vcf", "_summary.html"))

## ========
# EXAMPLE INPUT/ OUTPUT 
## ========
'''
User input should look something like this w/ argparse arguments
'''
#bash
##python compare_fasta_snps.py \
  ###--ref-dir user-data/Arabidopsis_thaliana/ncbi_dataset/data \
  ###--query-dir user-data/Arabidopsis_lyrata/ncbi_dataset/data \
  ###--out genomes/athal_vs_alyr_by_order.vcf \
  ###--mode order \
  ###--top-n 5








