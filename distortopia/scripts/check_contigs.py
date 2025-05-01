import os
import pandas as pd

# Get absolute path to this script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Correctly resolve paths relative to script location
snp_path = os.path.join(script_dir, "..", "genomes", "par_alignments", "snp_positions.tsv")
thaliana_path = os.path.join(script_dir, "..", "genomes", "par_alignments", "thaliana_contigs.txt")
lyrata_path = os.path.join(script_dir, "..", "genomes", "par_alignments", "lyrata_contigs.txt")

# Load data
tsv = pd.read_csv(snp_path, sep="\t")
ref_ids = set(open(thaliana_path).read().splitlines())
alt_ids = set(open(lyrata_path).read().splitlines())

# Identify missing
missing = []
for ref, alt in zip(tsv["Query"], tsv["Target"]):
    if ref not in ref_ids or alt not in alt_ids:
        missing.append((ref, alt))

# Print results
print(f"Total missing contig pairs: {len(missing)}")
for pair in missing[:10]:
    print("Missing:", pair)