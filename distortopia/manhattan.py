import matplotlib.pyplot as plt
import pandas as pd

'''
This script reads a .paf file, counts SNPs from the `cs` tag,
bins them by 10kb windows, and plots all scaffolds in a single
pseudochromosome layout. Scaffolds are concatenated along the x-axis.
Only high-SNP scaffolds are highlighted in red.
'''

# === PARSE PAF AND BIN SNPS ===
def parse_paf_for_snps(paf_path, bin_size=10000):
    snp_data = []

    with open(paf_path) as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 12:
                continue

            target = cols[5]           # scaffold name
            t_start = int(cols[7])    # start position on target
            cs_tag = [c[5:] for c in cols[12:] if c.startswith("cs:Z:")]

            if cs_tag:
                cs = cs_tag[0]
                snps = cs.count("*")  # SNPs are marked as substitutions
                bin_pos = (t_start // bin_size) * bin_size
                snp_data.append((target, bin_pos, snps))

    df = pd.DataFrame(snp_data, columns=["scaffold", "bin", "snps"])
    binned = df.groupby(["scaffold", "bin"])["snps"].sum().reset_index()
    return binned

# === PLOT PSEUDOCHROMOSOME LAYOUT ===
def plot_pseudochromosome(binned_df, highlight_threshold=250):
    offset = 0
    all_data = []
    scaffold_offsets = {}
    colors = []

    # choose colors
    color_normal = 'lightgray'
    color_highlight = 'red'

    # shift bins for each scaffold to get continuous x-axis
    for scaffold in binned_df["scaffold"].unique():
        sub = binned_df[binned_df["scaffold"] == scaffold].copy()
        sub["x"] = sub["bin"] + offset
        all_data.append(sub)
        scaffold_offsets[scaffold] = offset
        max_snp = sub["snps"].max()

        # red if max SNPs in a bin is high
        if max_snp >= highlight_threshold:
            colors.append(color_highlight)
        else:
            colors.append(color_normal)

        offset += sub["bin"].max() + 10000  # add gap between scaffolds

    final = pd.concat(all_data)
    fig, ax = plt.subplots(figsize=(12, 4))

    for i, (scaffold, subdf) in enumerate(final.groupby("scaffold")):
        ax.bar(subdf["x"], subdf["snps"], width=10000,
               color=colors[i], edgecolor='black')

    ax.set_title("Genome-wide SNP Density (pseudochromosome layout)")
    ax.set_xlabel("Pseudochromosome Position (bp)")
    ax.set_ylabel("SNPs per 10kb")
    plt.tight_layout()
    plt.show()

# === MAIN ===
if __name__ == "__main__":
    paf_file = "genomes/athal_vs_alyr.paf"
    binned = parse_paf_for_snps(paf_file)
    plot_pseudochromosome(binned)
