import os
import glob
import subprocess
from Bio import SeqIO
from io import StringIO #to parse VCF text input into pandas 
import argparse #to establish command-line flags for users
import pandas as pd #to genrate HTML with color-codes SNPs

"""
This script allows a user to select a genome FASTA from each of two species,
then compares them base-by-base for all matching contigs (by order and/or length).

It outputs a real, biologically meaningful VCF of SNPs between the two genome assemblies.
"""
def choose_fasta(species_name, species_dir):
    fasta_paths = glob.glob(os.path.join(species_dir, "ncbi_dataset", "data", "*", "*.fna*"), recursive=True)
    if not fasta_paths:
        raise FileNotFoundError(f"No FASTA files found for {species_name} in {species_dir}") #if there are no FASTA files found from NCBI

    print(f"\n[SELECT FASTA FOR {species_name}]") #giving user option to select which species
    for i, path in enumerate(fasta_paths):
        print(f"{i+1}. {path}")

    while True:
        try:
            index = int(input(f"Enter number [1-{len(fasta_paths)}]: ")) - 1 #user inputs number from 1-(some num of species fasta files)
            #NCBI might have multiple .fna files, hence the selection option 
            if 0 <= index < len(fasta_paths):
                return fasta_paths[index]
            else:
                print("Invalid selection. Try again.")
        except ValueError:
            print("Please enter a valid number.")

def run_minimap2(ref_fasta, query_fasta, output_vcf):
    script_dir = os.path.dirname(__file__)
    minimap_path = "minimap2"
    paf_path = output_vcf.replace(".vcf", ".paf")

    print("\nRunning minimap2...")
    with open(paf_path, "w") as paf_out:
        subprocess.run(
            [minimap_path, "-cx", "asm5", ref_fasta, query_fasta],
            stdout=paf_out,
            check=True
        )

    ref_seqs = {rec.id: str(rec.seq) for rec in SeqIO.parse(ref_fasta, "fasta")}

    print("Calling variants from PAF...")
    with open(paf_path) as paf, open(output_vcf, "w") as vcf:
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write("##source=minimap2+python\n")
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for line in paf:
            cols = line.strip().split("\t")
            if len(cols) < 12 or not cols[0] or not cols[5]:
                continue

            qname, qstart, qend = cols[0], int(cols[2]), int(cols[3])
            tname, tstart, tend = cols[5], int(cols[7]), int(cols[8])
            ref_seq = ref_seqs.get(tname)
            if not ref_seq:
                continue
            ref_sub = ref_seq[tstart:tend]

            for field in cols[12:]:
                if field.startswith("cs:Z:"):
                    cs = field[5:]
                    ref_pos = tstart
                    i = 0
                    while i < len(cs):
                        if cs[i] == ":":
                            i += 1
                            num = ""
                            while i < len(cs) and cs[i].isdigit():
                                num += cs[i]
                                i += 1
                            ref_pos += int(num)
                        elif cs[i] == "*":
                            if i + 2 < len(cs):
                                ref_base = cs[i+1]
                                alt_base = cs[i+2]
                                vcf.write(f"{tname}\t{ref_pos+1}\t.\t{ref_base}\t{alt_base}\t.\tPASS\t.\n")
                                ref_pos += 1
                            i += 3
                        else:
                            i += 1

def gen_HTML(vcf_path, html_out = "snp_summary.html"):
    """
    Loads a user-generated VCF file, summarizes SNP per contig, exports a styled HTML table.
    """ 
    with open(vcf_path) as f:
        lines = [line for line in f if not line.startswith("#")]

    df = pd.read_csv(
        StringIO("".join(lines)),
        sep="\t", #separating columns 
        names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"] #column names
    )

    # Filter biallelic SNPs
    valid = {"A", "C", "G", "T"} #looking for base pairs
    df = df[(df["REF"].isin(valid)) & (df["ALT"].isin(valid)) & (df["REF"] != df["ALT"])]

    # Count SNPs per contig
    summary = df["CHROM"].value_counts().reset_index()
    summary.columns = ["CHROM", "SNP_COUNT"]

    # Save to HTML with simple formatting
    summary.to_html(html_out, index=False)
    print(f"SNP summary table written to: {html_out}")

def parse_args():
    parser = argparse.ArgumentParser(description="Compare top contigs from two genomes by length or order")
    parser.add_argument("--ref-dir", required=True, help="Directory containing reference FASTAs")
    parser.add_argument("--query-dir", required=True, help="Directory containing query FASTAs")
    parser.add_argument("--out", required=True, help="Output VCF file path")
    parser.add_argument("--mode", choices=["length", "order"], default="length",
                        help="How to match contigs: by length (default) or order")
    parser.add_argument("--top-n", type=int, default=5, help="Number of contigs to compare (default: 5)")
    parser.add_argument("--summary", action="store_true", help="Generate an HTML summary table of SNP counts per contig")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    ref_fasta = choose_fasta("Species 1", args.ref_dir)
    query_fasta = choose_fasta("Species 2", args.query_dir)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    run_minimap2(ref_fasta, query_fasta, args.out)

    if args.summary:
        html_out = args.out.replace(".vcf", "_summary.html")
        gen_HTML(args.out, html_out)

## ========
# EXAMPLE INPUT/OUTPUT 
## ========
'''
User input should look something like this w/ argparse flags:
'''
#bash
##python compare_fasta_snps.py \
  ###--ref-dir user-data/Arabidopsis_thaliana/ncbi_dataset/data \
  ###--query-dir user-data/Arabidopsis_lyrata/ncbi_dataset/data \
  ###--out genomes/athal_vs_alyr_by_order.vcf \
  ###--mode order \
  ###--top-n 5

'''
User output should be a VCF and an HTML (optional). 
To generate the HTML user will need to use the --summary flag. 
To view/open the HTML user should enter this command in their terminal:
'''
#bash
##open genomes/athal_vs_alyr_by_order_summary.html





