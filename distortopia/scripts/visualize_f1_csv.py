import pandas as pd #for loading and manipulating the recombination table
import plotly.graph_objects as go #for interactive plotting with HTML output

#Load the F1 recombination table
df = pd.read_csv("genomes/f1_simulations/f1_table.csv")

#Convert 'crossover_pos' to numeric values
#turning "None" into NaN for easy logic
df["crossover_pos"] = pd.to_numeric(df["crossover_pos"], errors="coerce")

#Define colors for each parental origin
color_map = {"ref": "royalblue", "alt": "orangered"}

#Create a unique label for each row combining contig ID and replicate number
df["label"] = df["ref_chrom"] + "_rep" + df["rep"].astype(str)

#Sort table so that contigs are grouped together and replicates appear in order
df.sort_values(by=["ref_chrom", "rep"], inplace=True)

#Initialize blank interactive plot
fig = go.Figure()

#Track legend entries to avoid duplicating them
seen_legends = set()

#Loop through each row in the recombination table
for i, row in enumerate(df.itertuples()):
    y = row.label  #Y-axis label = contig + replicate
    #If crossover is missing, assign an arbitrary long segment (used to show single-color inheritance)
    crossover = row.crossover_pos if pd.notna(row.crossover_pos) else 1e7

    # ----- LEFT BLOCK -----
    # Only show the legend once per parental origin
    show_left = row.left not in seen_legends
    # Add a horizontal bar for the left segment of the chromosome
    fig.add_trace(go.Bar(
        x=[crossover],                #Width of the bar
        y=[y],                        #Y position
        name=row.left,               #Legend label
        marker_color=color_map[row.left],  #Color from map
        orientation='h',             #Horizontal bars
        hovertemplate=f"{y}<br>Left: {row.left}<br>Length: {int(crossover):,}<extra></extra>",
        showlegend=show_left         #Only show legend if not shown before
    ))
    seen_legends.add(row.left)       #Mark this parent as "seen"

    # ----- RIGHT BLOCK (if crossover occurred) -----
    if pd.notna(row.crossover_pos):
        show_right = row.right not in seen_legends
        # Add a horizontal bar for the right segment
        fig.add_trace(go.Bar(
            x=[1e7],#Arbitrary long segment (to extend past crossover)
            y=[y],
            base=[crossover], #Start at the crossover position
            name=row.right,
            marker_color=color_map[row.right],
            orientation='h',
            hovertemplate=f"{y}<br>Right: {row.right}<br>Start: {int(crossover):,}<extra></extra>",
            showlegend=show_right
        ))
        seen_legends.add(row.right)

# ----- Plot Layout -----
fig.update_layout(
    barmode='stack',  #Stack left and right segments on the same row
    title="F1 Recombination Map: All Replicates (Grouped by Contig)",
    xaxis_title="Position on Contig",
    yaxis_title="Contig ID (Replicate)",
    yaxis=dict(automargin=True),    #Prevent y-labels from getting clipped
    legend_title="Parental Origin",
    template="plotly_white",        #Clean plot
    height=20 * len(df) + 200,      #Dynamic plot height based on number of rows
    width=1000                      #Fixed plot width
)

#Save to an HTML file
fig.write_html("f1_recombination_dashboard.html")

# Notify user
print("Dashboard written to f1_recombination_dashboard.html")
