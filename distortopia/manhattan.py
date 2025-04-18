import matplotlib.pyplot as plt
import pandas as pd
import re
'''
This script generates a genome-wide variant density/positional map. 
It uses the .paf file, and associates SNP positions and contigs. 
It requires users have: (1) .paf file containing target contig names and (2) SNP positions 
or approximate locations. These data will allow users to plot each variant along the 
genomic x-axis, grouped by contig.

The script groups SNPs by scaffold (e.g. chromosome/contig identifiers from .paf file) and
bins them into 10kb windows for plotting SNP density. Each scaffold is plotted separately 
for clarity.
'''

def parse_paf_for_snps(paf_path, bin_size=10000):
	"""
	Parses user-provided .paf file and counts SNPs into 10kb bins
	"""
    data = []

    with open(paf_path) as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 12:
                continue

            query = cols[0]
            target = cols[5]
            target_start = int(cols[7])
            target_end = int(cols[8])
            aln_len = int(cols[11])

            #Extract cs tags (variant annotation)
            cs_tag = ""
            for col in cols[12:]:
                if col.startswith("cs:Z:"):
                    cs_tag = col[5:]
                    break

            #Count SNPs using cs tag
            snp_count = cs_tag.count("*")

            #Store for every alignment
            data.append({
                "scaffold": target,
                "start": target_start,
                "end": target_end,
                "snps": snp_count
            })

    df = pd.DataFrame(data)

    # dd bin column (bin start coordinate)
    df["bin"] = (df["start"] // bin_size) * bin_size

    #Group by scaffold and bin, sum SNPs
    binned = df.groupby(["scaffold", "bin"])["snps"].sum().reset_index()
    return binned

def plot_snp_density(binned_df):
    scaffolds = binned_df["scaffold"].unique()

    for scaffold in scaffolds:
        sub = binned_df[binned_df["scaffold"] == scaffold]

        plt.figure(figsize=(10, 4))
        plt.bar(sub["bin"], sub["snps"], width=10000, color='teal', edgecolor='black')
        plt.title(f"SNP Density - {scaffold}")
        plt.xlabel("Position (bp)")
        plt.ylabel("SNP count (per 10kb)")
        plt.tight_layout()
        plt.show()

# === MAIN EXECUTION ===
if __name__ == "__main__":
	
    paf_file = "genomes/athal_vs_alyr.paf"  # <- Adjust path if needed
    binned_snps = parse_paf_for_snps(paf_file)
    plot_snp_density(binned_snps)


