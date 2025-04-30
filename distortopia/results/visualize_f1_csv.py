import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV
df = pd.read_csv("genomes/f1_simulations/f1_table.csv")

# Set up the color map
color_map = {'ref': '#4C72B0', 'alt': '#DD8452'}

#Group by F1 replicate
for rep, group in df.groupby("rep"):
    fig, ax = plt.subplots(figsize=(10, 5))
    y_labels = []

    for i, row in enumerate(group.itertuples()):
        # Use contig index for vertical placement
        y = len(group) - i
        contig = row.ref_chrom
        y_labels.append(contig)

        # Plot left block
        ax.broken_barh([(0, row.crossover_pos if row.crossover_pos != "None" else 1e7)],
                       (y - 0.3, 0.6),
                       facecolors=color_map[row.left])

        #Plot right block (if crossover happened)
        if row.crossover_pos != "None":
            ax.broken_barh([(int(row.crossover_pos), 1e7)],
                           (y - 0.3, 0.6),
                           facecolors=color_map[row.right])

            #Draw vertical line at crossover position
            ax.axvline(int(row.crossover_pos), ymin=0.1, ymax=0.9, color="black", linestyle="--", lw=1)

    ax.set_yticks(range(1, len(group) + 1))
    ax.set_yticklabels(y_labels[::-1])
    ax.set_xlabel("Position on Chromosome")
    ax.set_title(f"F1 Hybrid Recombination Map – Replicate {rep}")
    ax.set_xlim(0, 2e7)  # adjust as needed based on your genome
    plt.tight_layout()
    plt.savefig(f"f1_rep{rep}_recombination_plot.png")
    plt.show()