import os
import glob
import subprocess
import time
from Bio import SeqIO #handles biological sequence files (not really used here)
import argparse #for command-line argument parsing
import pandas as pd #for data analysis and manipulation

def choose_fasta(species_name, species_dir):
    """
    Lists all FASTA files in the species NCBI directory and prompts user to select.
    """
    fasta_paths = glob.glob(os.path.join(species_dir, "**", "*.fna*"), recursive=True) #finds/fetches all .fna files 
    if not fasta_paths: #raising error!!
        raise FileNotFoundError(f"No FASTA files found for {species_name} in {species_dir}")
    print(f"\n[SELECT FASTA FOR {species_name}]")
    for i, path in enumerate(fasta_paths):
        print(f"{i+1}. {path}") #displays all available FASTA files
    #prompts user to sleect a FASTA file from list 
    while True:
        try:
            index = int(input(f"Enter number [1-{len(fasta_paths)}]: ")) - 1
            if 0 <= index < len(fasta_paths):
                return fasta_paths[index]
            else:
                print("Invalid selection. Try again.") #if user tries to enter invalid index
        except ValueError:
            print("Please enter a valid number.")

def run_minimap2(ref_fasta, alt_fasta, paf_path, threads=4, preset="asm5"):
    """
    Run minimap2 to align two genome FASTA files and write output to .paf file.
    """
    minimap_path = "/opt/homebrew/bin/minimap2" #path to minimap2 binary file
    print(f"Running minimap2 with {threads} threads...")
    start_time = time.time() #start timer for timing alignment process
    #run minimap2 and write output to a .paf file
    with open(paf_path, "w") as paf_out:
        subprocess.run(
            [minimap_path, "-t", str(threads), "-cx", preset, "--cs=short", ref_fasta, alt_fasta],
            stdout=paf_out,
            check=True
        )
    elapsed_time = time.time() - start_time
    print(f"Finished minimap2 in {elapsed_time:.2f} seconds.")

def summary_of_paf(paf_path, html_out="alignment_summary.html", snp_out="snp_positions.tsv"):
    """
    Parse the .paf file and extract alignment stats + SNP positions.
    """
    results = [] #store parsed alignment data
    with open(paf_path) as f:
        for line in f:
            cols = line.strip().split("\t") #splitting lines 
            if len(cols) < 12:
                continue #skip any incomplete lines 
            #basic alignment info 
            query, target = cols[0], cols[5]
            match_len = int(cols[10]) #number of residue matches
            aln_len = int(cols[11]) #alignment block length
            #extract cs:Z tag for *, -, +
            cs_tag = ""
            for col in cols[12:]:
                if col.startswith("cs:Z:"):
                    cs_tag = col[5:]
                    break
            #counters for counting alignment features 
            matches, snps, indels = 0, 0, 0 #initializing counters
            ref_pos = int(cols[7]) #start position on reference genome
            snp_positions = [] #initializing list of SNP positions
            #parses cs:Z tags 
            i = 0 #initializing index 
            while i < len(cs_tag):
                if cs_tag[i] == ":": #exact-bair pair matches
                    i += 1
                    num = ""
                    while i < len(cs_tag) and cs_tag[i].isdigit():
                        num += cs_tag[i]
                        i += 1
                    matches += int(num)
                    ref_pos += int(num)
                elif cs_tag[i] == "*": #SNPs
                    snps += 1
                    snp_positions.append(ref_pos)
                    ref_pos += 1
                    i += 3 #skip over any mismatched base-pairs
                elif cs_tag[i] == "+": #insertions
                    indels += 1
                    i += 1
                    while i < len(cs_tag) and cs_tag[i].isalpha():
                        i += 1
                elif cs_tag[i] == "-": #deletions
                    indels += 1
                    i += 1
                    del_len = 0
                    while i < len(cs_tag) and cs_tag[i].isalpha():
                        del_len += 1
                        i += 1
                    ref_pos += del_len
                else:
                    i += 1
            #store results for each alignment 
            results.append({ #adding column headers
                "Query": query, 
                "Target": target,
                "Map_q": aln_len,
                "Matches": matches,
                "SNPs": snps,
                "Indels": indels,
                "SNP_Positions": ", ".join(map(str, snp_positions))
            })
    #convert results to a df and then save to a .tsv
    df = pd.DataFrame(results)
    df[["Query", "Target", "SNP_Positions"]].to_csv(snp_out, sep="\t", index=False)

    #styling dataframe for HTML output summary table
    styled = df.style \
    .background_gradient(subset=["SNPs", "Indels"], cmap="Reds") \
    .highlight_max(color="lightgreen", axis=0, subset=["Matches"]) \
    .set_caption("Minimap2 Alignment Summary") \
    .set_table_styles([
        {'selector': 'th', 'props': [('background-color', '#f2f2f2'), ('color', '#333'), ('font-size', '12px')]},
        {'selector': 'caption', 'props': [('caption-side', 'top'), ('font-size', '16px'), ('font-weight', 'bold')]}
    ])

#Write the HTML file to disk
    with open(html_out, "w") as f:
        f.write(styled.to_html())

    print(f"HTML summary written to: {html_out}")
    print(f"SNP table written to: {snp_out}")

def parse_args():
    parser = argparse.ArgumentParser(description="Align two genomes and extract SNP summary")
    parser.add_argument("--ref-dir", required=True, help="Path to reference genome directory")
    parser.add_argument("--alt-dir", required=True, help="Path to alternate genome directory")
    parser.add_argument("--outdir", default="genomes/alignment", help="Directory to store alignment outputs")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads for minimap2")
    parser.add_argument("--preset", default="asm5", help="Minimap2 preset (e.g., asm5, asm10)")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    #prompts user to choose a ref and alt FASTA file
    ref_fasta = choose_fasta("Parent 1", args.ref_dir)
    alt_fasta = choose_fasta("Parent 2", args.alt_dir)
    #creates an output directory folder if one doesn't already exist 
    os.makedirs(args.outdir, exist_ok=True)

    #defining file output paths 
    paf_path = os.path.join(args.outdir, "ref_vs_alt.paf")
    html_out = os.path.join(args.outdir, "ref_vs_alt_summary.html")
    snp_out = os.path.join(args.outdir, "snp_positions.tsv")

    run_minimap2(ref_fasta, alt_fasta, paf_path, threads=args.threads, preset=args.preset)
    summary_of_paf(paf_path, html_out, snp_out)

    #automatically opens HTML (default web browser)
    os.system(f"open {html_out}")

#=====
#EXAMPLE CLI
#=====

#(base) riyarampalli@Riyas-MacBook-Pro-7 distortopia % python scripts/align_parents.py \
  #--ref-dir input-data/Arabidopsis_thaliana \
  #--alt-dir input-data/Arabidopsis_lyrata \
  #--outdir genomes/par_alignments \
  #--threads 8




