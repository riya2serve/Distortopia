import argparse
import os
import subprocess
import pandas as pd

def run_minimap2(ref_fasta, reads_fasta, paf_out, threads=4):
    cmd = [
        "minimap2", "-x", "map-pb", "-t", str(threads),
        ref_fasta, reads_fasta
    ]
    with open(paf_out, "w") as out_f:
        subprocess.run(cmd, stdout=out_f, check=True)

def parse_paf(paf_path):
    cols = [
        "ref_name", "ref_len", "ref_start", "ref_end", "strand",
        "alt_name", "alt_len", "alt_start", "alt_end",
        "n_match", "aln_len", "mapq"
    ]
    records = []
    with open(paf_path) as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) >= 12:
                record = {col: fields[i] for i, col in enumerate(cols)}
                records.append(record)
    return pd.DataFrame(records)

def main(ref_fasta, alt_fasta, snp_tsv, output_dir, reps, threads):
    os.makedirs(output_dir, exist_ok=True)

    for rep in range(1, reps + 1):
        reads_path = os.path.join(output_dir, f"f1_reads_rep{rep}.fasta")

        for label, parent_fasta in [("ref1", ref_fasta), ("ref2", alt_fasta)]:
            paf_out = os.path.join(output_dir, f"rep{rep}_vs_{label}.paf")
            print(f"Mapping rep{rep} reads to {label}...")
            run_minimap2(parent_fasta, reads_path, paf_out, threads)

            # Parse PAF and generate HTML summary
            df = parse_paf(paf_out)
            if not df.empty:
                df["aln_len"] = pd.to_numeric(df["aln_len"])
                df["mapq"] = pd.to_numeric(df["mapq"])
                html_path = os.path.join(output_dir, f"rep{rep}_vs_{label}_summary.html")
                df[["query_name", "target_name", "aln_len", "mapq"]].to_html(html_path, index=False)
                print(f"Saved HTML summary: {html_path}")
            else:
                print(f"No alignments found in {paf_out}")

def parse_args():
    parser = argparse.ArgumentParser(description="Map simulated long reads to parents and summarize alignments.")
    parser.add_argument("--ref-fasta", required=True, help="Parent 1 FASTA")
    parser.add_argument("--alt-fasta", required=True, help="Parent 2 FASTA")
    parser.add_argument("--snp-tsv", required=True, help="SNP TSV (not used here)")
    parser.add_argument("--outdir", default="genomes/f1_simulations", help="Where long reads and PAFs live")
    parser.add_argument("--reps", type=int, default=1, help="Number of replicates")
    parser.add_argument("--threads", type=int, default=4, help="Threads for minimap2")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    main(args.ref_fasta, args.alt_fasta, args.snp_tsv, args.outdir, args.reps, args.threads)
