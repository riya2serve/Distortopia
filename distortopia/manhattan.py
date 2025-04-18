import matplotlib.pyplot as plt
import pandas as pd

'''
This script generates a genome-wide SNP density map from a .paf file.
It bins SNPs into 10kb windows based on their position on each scaffold (contig).
The result is one density plot per scaffold, showing how SNPs are distributed.
'''

import matplotlib.pyplot as plt
import pandas as pd

'''
This script generates a genome-wide SNP density map from a .paf file.
It bins SNPs into 10kb windows based on their position on each scaffold (contig).
The result is one density plot per scaffold, showing how SNPs are distributed.
'''

# === PARSE .PAF FILE AND COUNT SNPs PER BIN ===
def parse_paf_for_snps(paf_path, bin_size=10000):
    snp_data = []

    with open(paf_path) as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 12:
                continue

            target = cols[5]           # scaffold name
            t_start = int(cols[7])    # alignment start position on target
            cs_tag = [c[5:] for c in cols[12:] if c.startswith("cs:Z:")]

            if cs_tag:
                cs = cs_tag[0]
                snps = cs.count("*")  # count SNPs as substitutions (*)

                bin_pos = (t_start // bin_size) * bin_size
                snp_data.append((target, bin_pos, snps))

    # convert to DataFrame
    df = pd.DataFrame(snp_data, columns=["scaffold", "bin", "snps"])
    binned = df.groupby(["scaffold", "bin"])["snps"].sum().reset_index()
    return binned

# === PLOT SNP DENSITY PER SCAFFOLD ===
def plot_snp_density(binned_df):
    scaffolds = binned_df["scaffold"].unique()

    for scaffold in scaffolds:
        sub = binned_df[binned_df["scaffold"] == scaffold]

        plt.figure(figsize=(10, 4))
        plt.bar(sub["bin"], sub["snps"], width=10000, color='teal', edgecolor='black')
        plt.title(f"SNP Density - {scaffold}")
        plt.xlabel("Position (bp)")
        plt.ylabel("SNPs per 10kb")
        plt.tight_layout()
        plt.show()

# === MAIN ===
if __name__ == "__main__":
    paf_file = "genomes/athal_vs_alyr.paf"  # <- update path if needed
    binned = parse_paf_for_snps(paf_file)
    plot_snp_density(binned)


