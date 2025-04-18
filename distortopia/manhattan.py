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

    # pick colors
    color_normal = 'gray'
    color_highlight = 'red'

    # Sort scaffolds by SNP count (helps avoid clutter)
    scaffold_order = binned_df.groupby("scaffold")["snps"].sum().sort_values(ascending=False).index

    fig, ax = plt.subplots(figsize=(14, 5))  # Wider layout

    for i, scaffold in enumerate(scaffold_order):
        sub = binned_df[binned_df["scaffold"] == scaffold].copy()
        sub["x"] = sub["bin"] + offset
        all_data.append(sub)

        max_snp = sub["snps"].max()
        color = color_highlight if max_snp > highlight_threshold else color_normal

        # plot bar
        ax.bar(sub["x"], sub["snps"], width=10000, color=color, edgecolor='black')

        # add vertical line to separate scaffolds
        offset += sub["bin"].max() + 50000  # more spacing between scaffolds

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
