import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Load the generated CSV from genomes folder
df = pd.read_csv("genomes/f1_simulations/f1_table.csv")
df["crossover_pos"] = pd.to_numeric(df["crossover_pos"], errors="coerce")  # turns "None" into NaN

# Set up the color map
color_map = {'ref': '#4C72B0', 'alt': '#DD8452'}

# Group by F1 replicate
for rep, group in df.groupby("rep"):
    group = group.drop_duplicates(subset=["rep", "ref_chrom"])
    fig, ax = plt.subplots(figsize=(10, 5))
    y_labels = []

    for i, row in enumerate(group.itertuples()):
        y = len(group) - i
        contig = row.ref_chrom
        y_labels.append(contig)

        crossover = row.crossover_pos if not pd.isna(row.crossover_pos) else 1e7

        # Plot left block
        ax.broken_barh([(0, crossover)],
                       (y - 0.3, 0.6),
                       facecolors=color_map[row.left])

        # Plot right block if crossover occurred
        if not pd.isna(row.crossover_pos):
            ax.broken_barh([(row.crossover_pos, 1e7)],
                           (y - 0.3, 0.6),
                           facecolors=color_map[row.right])
            ax.axvline(row.crossover_pos, ymin=0.1, ymax=0.9,
                       color="black", linestyle="--", lw=1)

    ax.set_yticks(range(1, len(group) + 1))
    ax.set_yticklabels(y_labels[::-1])
    ax.set_xlabel("Position on Chromosome")
    ax.set_title(f"F1 Hybrid Recombination Map – Replicate {rep}")
    ax.set_xlim(0, 2e7)

    # Add legend
    legend_elements = [
        Patch(facecolor=color_map['ref'], label='ref'),
        Patch(facecolor=color_map['alt'], label='alt')
    ]
    ax.legend(handles=legend_elements, loc='lower right')

    plt.tight_layout()
    plt.savefig(f"f1_rep{rep}_recombination_plot.png")
    plt.show()
