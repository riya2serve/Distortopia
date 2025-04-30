import pandas as pd #to create dataframe
import plotly.graph_objects as go #to generate interactive HTML

#Load CSV
df = pd.read_csv("genomes/f1_simulations/f1_table.csv")
df["crossover_pos"] = pd.to_numeric(df["crossover_pos"], errors="coerce")

# Color map
color_map = {"ref": "royalblue", "alt": "orangered"}

#Plot one interactive figure per replicate
for rep, group in df.groupby("rep"):
    group = group.drop_duplicates(subset=["rep", "ref_chrom"])
    fig = go.Figure()
    y_labels = []
    max_x = 0

    for i, row in enumerate(group.itertuples()):
        y = len(group) - i
        y_labels.append(row.ref_chrom)
        crossover = row.crossover_pos if pd.notna(row.crossover_pos) else 1e7
        max_x = max(max_x, crossover if pd.notna(row.crossover_pos) else 1e7)

        #Left block
        fig.add_trace(go.Bar(
            x=[crossover],
            y=[row.ref_chrom],
            name=f"{row.left}",
            marker_color=color_map[row.left],
            orientation='h',
            hovertemplate=f"{row.ref_chrom}<br>Left: {row.left}<br>Length: {int(crossover):,}<extra></extra>"
        ))

        #Right block
        if pd.notna(row.crossover_pos):
            fig.add_trace(go.Bar(
                x=[1e7],  # long enough to fill right side
                y=[row.ref_chrom],
                base=[crossover],
                name=f"{row.right}",
                marker_color=color_map[row.right],
                orientation='h',
                hovertemplate=f"{row.ref_chrom}<br>Right: {row.right}<br>Start: {int(crossover):,}<extra></extra>"
            ))

    fig.update_layout(
        barmode='stack',
        title=f"F1 Hybrid Recombination Map – Replicate {rep}",
        xaxis_title="Position on Contig", #xaxis title
        yaxis_title="Contig ID", #yaxis title
        yaxis=dict(categoryorder="array", categoryarray=list(reversed(y_labels))),
        legend_title="Parental Origin",
        template="plotly_white",
        height=30 * len(group) + 100,
        width=1000
    )

    # Save as interactive HTML
    fig.write_html(f"f1_rep{rep}_recombination_plot.html") #naming file 
    print(f"Saved: f1_rep{rep}_recombination_plot.html")
