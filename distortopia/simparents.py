import os #interacts with user's operating system 
import glob #used for finding files/patterns
import subprocess #allows users to run shell commands 
import time #to track run time
from Bio import SeqIO #BioPython module for reading and writing sequence files 
import argparse #allows for command-line arguments/flags from/for users
import pandas as pd #for working with dataframes and generated HTML style summaries

"""
This script allows a user to select a genome FASTA from each of two species, then compares them. 
It outputs a biologically meaningful summary of SNPs, indels, and alignment metrics in an HTML summary table.
"""
def choose_fasta(species_name, species_dir):
    """
    Lists all FASTA files in the species NCBI directory and prompts user to select.
    """
    #finds all fasta files ending with .fna (or similar) inside the nested 'ncbi_dataset' directory 
    fasta_paths = glob.glob(os.path.join(species_dir, "ncbi_dataset", "data", "*", "*.fna*"), recursive=True)
    if not fasta_paths: #raises an error if no FASTA files are found 
        raise FileNotFoundError(f"No FASTA files found for {species_name} in {species_dir}")  
    #displays list of all FASTA files avaiable for user to select from 
    print(f"\n[SELECT FASTA FOR {species_name}]") 
    for i, path in enumerate(fasta_paths):
        print(f"{i+1}. {path}") #displays numbers (1 - n); user selects from these to choose FASTA file 

    while True: #loop -- user must select a valid index between (1-n)
        try:
            index = int(input(f"Enter number [1-{len(fasta_paths)}]: ")) - 1 #user input 
            if 0 <= index < len(fasta_paths):
                return fasta_paths[index] #returns selected FASTA file path 
            else:
                print("Invalid selection. Try again.") #if user selects out-of-range index number 
        except ValueError:
            print("Please enter a valid number.") #handles non-integer input from users 

def run_minimap2(ref_fasta, query_fasta, paf_path, threads = 4, preset = "asm5"):
    """
    Runs minimpa2 to align two genome FASTA files and writes output to .paf file
    """
    minimap_path = "/opt/homebrew/bin/minimap2" #the full path to the minimap2 binary source file 

    print(f"Running minimap2 with {threads} threads......") #prints message to let user know minimap2 is running
    start_time = time.time() #starts timer (allow user to test how long minimap2 takes to run)

    #run minimap2 with specified user settings and write output file:
    with open(paf_path, "w") as paf_out:
        subprocess.run(
            [minimap_path, "-t", str(threads), "-cx", preset, "--cs=short", ref_fasta, query_fasta], #user-command
            stdout=paf_out, #redirects minimap2 to user-selected output file 
            check=True #raises an exception if minimap2 fails to run 
        )
    #calculates how long minimap2 took to run
    elapsed_time = time.time() - start_time #end time - start time 
    print(f"Finished minimap2 in {elapsed_time:.2f} seconds.. Now loading reference sequences...")

def summary_of_paf(paf_path, html_out = "alignment_summary.html"):
    """
    Parses a .paf file to extract alignment metrics and variant types, and outputs a stylized HTML table.
    """
    results = [] #initializes list to collect/store parsed results

    #opens and reads .paf alignment file 
    with open(paf_path) as f:
        for line in f:
            cols = line.strip().split("\t") #splitting .paf file into tab separated fields 
            if len(cols) < 12:
                continue #skipping lines that don't contain expected fields 

            query, target = cols[0], cols[5] #extracts query_sequence and target_sequence names 
            match_len = int(cols[10]) #extracts num. of bases that matched between the two 
            aln_len = int(cols[11]) #extracts total alignment length 

            #searches for cs tag -- encodes for SNPS, indels, matches, etc. 
            cs_tag = ""
            for col in cols[12:]:
                if col.startswith("cs:Z:"): #encodes the alignment at the base-pair level
                    cs_tag = col[5:] #extracts content associated with cs tag 
                    break

            matches, snps, indels = 0, 0, 0 #initializing counters for these 
            i = 0
            while i < len(cs_tag):
                if cs_tag[i] == ":": #exact base-pair matches
                    # Exact match of length
                    i += 1
                    num = ""
                    while i < len(cs_tag) and cs_tag[i].isdigit():
                        num += cs_tag[i]
                        i += 1
                    matches += int(num) #match sequence length 
                elif cs_tag[i] == "*": #base-pair substitution
                    snps += 1 #one substitution 
                    i += 3  #skip ref + alt bases
                elif cs_tag[i] in "+-": #base-pair insertion or deletion
                    indels += 1 #one insertion or deletion 
                    i += 1
                    while i < len(cs_tag) and cs_tag[i].isalpha():
                        i += 1 #skips inserted/deleted sequences 
                else:
                    i += 1 #catches all
            #appending all metrics as a dictionary -- key:value
            results.append({
                "Query": query, #contig name for reference genome (in this case A. thaliana)
                "Target": target, #contig name from the query genome (in this case A. lyrata)
                "Map_q": aln_len, #mapping quality of base pairs that aligned
                "Matches": matches, #number of exact matching base pairs 
                "SNPs": snps, #number of single-nucleotide differences (in this case substitutions)
                "Indels": indels #number of insertions and deletions 
            })

    #converts dictionary results to pandas DataFrame
    df = pd.DataFrame(results)

    #applies color styling to emphasize SNP/indel levels + highest match
    styled = df.style\
        .background_gradient(subset=["SNPs", "Indels"], cmap="Reds")\
        .highlight_max(color="lightgreen", axis=0, subset=["Matches"])\
        .format({"SNP_rate": "{:.2%}", "Indel_rate": "{:.2%}"})\
        .set_caption("Minimap2 Alignment Summary")\
        .set_table_styles([
            {'selector': 'th', 'props': [('background-color', '#f2f2f2'), ('color', '#333'), ('font-size', '12px')]},
            {'selector': 'caption', 'props': [('caption-side', 'top'), ('font-size', '16px'), ('font-weight', 'bold')]}
        ])

    #Export styled summary table to HTML
    styled.to_html(html_out)
    print(f"Styled summary written to: {html_out}")

def parse_args():
    """
    Command-line argument interface.
    Users can provide paths, number of threads, and minimap2 settings.
    """
    parser = argparse.ArgumentParser(description="Align two genomes with minimap2 and summarize variant types")
    parser.add_argument("--ref-dir", required=True, help="Directory containing reference FASTAs")
    parser.add_argument("--query-dir", required=True, help="Directory containing query FASTAs")
    parser.add_argument("--out", required=True, help="Output prefix (used to create .paf and .html)")
    parser.add_argument("--threads", type=int, default=4, help="Number of CPUs to use with minimap2 (default: 4)")
    parser.add_argument("--preset", default = "asm5", help = "minimap2 preset (e.g., asm5, asm10, splice)")
    return parser.parse_args() #parse all arguments and return them 

if __name__ == "__main__":
    #Step 1. parse command-line arguments
    args = parse_args()

    #Step 2. prompt user to select one FASTA file for each species (interactively)
    ref_fasta = choose_fasta("Species 1", args.ref_dir)
    query_fasta = choose_fasta("Species 2", args.query_dir)
    os.makedirs(os.path.dirname(args.out), exist_ok=True) #makes output file if one doesn't exist already

    #Step 3: ensure output file directories exist 
    paf_path = args.out
    html_out = args.out.replace(".paf", "_summary.html")
    os.makedirs(os.path.dirname(paf_path), exist_ok=True) #makes output file if one doesn't exist already 

    #Step 4: Run minimap2 + summarize PAF
    run_minimap2(ref_fasta, query_fasta, paf_path, threads=args.threads, preset=args.preset)
    summary_of_paf(paf_path, html_out)

    #Step 5. open the HTML summary (automatically) in default browser
    os.system(f"open {html_out}")


## ==========
# EXAMPLE CLI 
## ==========
'''
User input should look something like this w/ argparse flags:
'''
#bash
##distortopia % python simparents.py \          
  ###--ref-dir user-data/Arabidopsis_thaliana \
  ###--query-dir user-data/Arabidopsis_lyrata \
  ###--out genomes/athal_vs_alyr.paf \
  ###--threads 8

"""
Recommend checking if your files are a reasonable size before committing to GitHub repo.
Run this command in your terminal (note: file path might be slighlty different for each user)
"""
#bash
##ls -lh genomes/athal_vs_alyr_summary.html

"""
Understand the cs:Z tag
"""
#bash
##cs:Z::50*at:12*gc-aaa+ttt:20
###: — exact match of length N
###*at — SNP (A in ref, T in query)
###-aaa — deletion in the query
###+ttt — insertion in the query

