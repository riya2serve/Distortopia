import pandas as pd

tsv = pd.read_csv("genomes/par_alignments/snp_positions.tsv", sep="\t")
ref_ids = set(open("thaliana_contigs.txt").read().splitlines())
alt_ids = set(open("lyrata_contigs.txt").read().splitlines())

missing_ref = [r for r in tsv["Query"] if r not in ref_ids]
missing_alt = [a for a in tsv["Target"] if a not in alt_ids]

print("Missing from thaliana:", len(missing_ref), "Examples:", missing_ref[:5])
print("Missing from lyrata:", len(missing_alt), "Examples:", missing_alt[:5])
