import os
import glob
import subprocess
import time #to track run time
import re 
from Bio import SeqIO
import argparse #to establish command-line flags for users
import pandas as pd #to genrate HTML with color-codes SNPs

"""
This script allows a user to select a genome FASTA from each of two species,
then compares them. It outputs a real, biologically meaningful summary table of assembly metrics.
"""
def choose_fasta(species_name, species_dir):
    """
    Lists all FASTA files in the species NCBI directory and prompts user to select.
    """
    fasta_paths = glob.glob(os.path.join(species_dir, "ncbi_dataset", "data", "*", "*.fna*"), recursive=True)
    if not fasta_paths:
        raise FileNotFoundError(f"No FASTA files found for {species_name} in {species_dir}")  # if there are no FASTA files found from NCBI

    print(f"\n[SELECT FASTA FOR {species_name}]")  # giving user option to select which species
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

def run_minimap2(ref_fasta, query_fasta, paf_path, threads = 4, preset = "asm5"):
    """
    Runs minimpa2 to align two genome FASTA files and writes output to .paf file
    """
    minimap_path = "/opt/homebrew/bin/minimap2"

    print(f"Running minimap2 with {threads} threads......")
    start_time = time.time()

    #run minimap2 and write output:
    with open(paf_path, "w") as paf_out:
        subprocess.run(
            [minimap_path, "-t", str(threads), "-cx", preset, "--cs=short", ref_fasta, query_fasta],
            stdout=paf_out,
            check=True
        )

    elapsed_time = time.time() - start_time
    print(f"Finished minimap2 in {elapsed_time:.2f} seconds.. Now loading reference sequences...")

def summary_of_paf(paf_path, html = "alignment_summary.html"):
    """
    Parses .paf file and generates a summary table of alignment metrics and variant types.
    """
    results = []

    with open(paf_path) as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 12:
                continue

            query, target = cols[0], cols[5]
            match_len = int(cols[10])
            aln_len = int(cols[11])

            # Look for cs tag to extract variant types
            cs_tag = ""
            for col in cols[12:]:
                if col.startswith("cs:Z:"):
                    cs_tag = col[5:]
                    break

            matches, snps, indels = 0, 0, 0
            i = 0
            while i < len(cs_tag):
                if cs_tag[i] == ":":
                    # Exact match of length
                    i += 1
                    num = ""
                    while i < len(cs_tag) and cs_tag[i].isdigit():
                        num += cs_tag[i]
                        i += 1
                    matches += int(num)
                elif cs_tag[i] == "*":
                    snps += 1
                    i += 3  # skip ref + alt bases
                elif cs_tag[i] in "+-":
                    indels += 1
                    i += 1
                    while i < len(cs_tag) and cs_tag[i].isalpha():
                        i += 1
                else:
                    i += 1

            results.append({
                "Query": query,
                "Target": target,
                "Aligned_bp": aln_len,
                "Matches": matches,
                "SNPs": snps,
                "Indels": indels
            })

    # Convert to pandas DataFrame and export as HTML
    df = pd.DataFrame(results)
    df.to_html(html_out, index=False)
    print(f"Summary written to: {html_out}")

def parse_args():
    """
    Command-line arguments
    """
    parser = argparse.ArgumentParser(description="Align two genomes with minimap2 and summarize variant types")
    parser.add_argument("--ref-dir", required=True, help="Directory containing reference FASTAs")
    parser.add_argument("--query-dir", required=True, help="Directory containing query FASTAs")
    parser.add_argument("--out", required=True, help="Output prefix (used to create .paf and .html)")
    parser.add_argument("--threads", type=int, default=4, help="Number of CPUs to use with minimap2 (default: 4)")
    parser.add_argument("--preset", default = "asm5", help = "minimap2 preset (e.g., asm5, asm10, splice)")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    #Step 1. Allow user to select FASTA files interactively
    ref_fasta = choose_fasta("Species 1", args.ref_dir)
    query_fasta = choose_fasta("Species 2", args.query_dir)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    #Step 2: Prepare file output paths
    paf_path = args.out.replace(".vcf", ".paf")
    html_out = args.out.replace(".vcf", "_summary.html")
    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    # Step 3: Run minimap2 + summarize PAF
    run_minimap2(ref_fasta, query_fasta, paf_path, threads=args.threads, preset=args.preset)
    summary_of_paf(paf_path, html_out)

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





