import pandas as pd
import argparse

# Extract mismatches from MD:Z and read seq
def parse_sam_to_vcf(sam_path, ref_label):
    records = []
    with open(sam_path, "r") as sam:
        for line in sam:
            if line.startswith("@"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 12:
                continue
            query, chrom, pos, mapq, seq, tags = fields[0], fields[2], int(fields[3]), fields[4], fields[9], fields[11:]
            md_tag = next((tag for tag in tags if tag.startswith("MD:Z:")), None)
            if not md_tag:
                continue
            md = md_tag.split(":")[-1]
            i = 0
            ref_pos = pos
            num = ""
            for char in md:
                if char.isdigit():
                    num += char
                else:
                    if num:
                        i += int(num)
                        ref_pos += int(num)
                        num = ""
                    if char.isalpha() and i < len(seq):
                        records.append({
                            "CHROM": chrom,
                            "POS": ref_pos,
                            "ID": query,
                            "REF": char,
                            "ALT": seq[i],
                            "QUAL": mapq,
                            "FILTER": "PASS",
                            "INFO": f"source={ref_label}"
                        })
                        i += 1
    return records

def main(sam_ref, sam_alt, output_path, rep):
    print(f"Parsing SAM files for rep {rep} to generate: {output_path}")
    combined = parse_sam_to_vcf(sam_ref, "ref") + parse_sam_to_vcf(sam_alt, "alt")
    pd.DataFrame(combined).to_csv(output_path, sep="\t", index=False)
    print("Done.")

def parse_args():
    parser = argparse.ArgumentParser(description="Generate VCF-like TSV from SAM alignments.")
    parser.add_argument("--sam-ref", required=True, help="SAM file of long reads mapped to reference genome")
    parser.add_argument("--sam-alt", required=True, help="SAM file of long reads mapped to alternative genome")
    parser.add_argument("--output", required=True, help="Output path for VCF-like TSV file")
    parser.add_argument("--rep", type=int, required=True, help="Replicate number")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    main(args.sam_ref, args.sam_alt, args.output, args.rep)


#========
# EXAMPLE CLI
#========

#python scripts/sam_to_vcf.py \
  #--sam-ref distortopia/genomes/f1_simulations/rep1_vs_ref.sam \
  #--sam-alt distortopia/genomes/f1_simulations/rep1_vs_alt.sam \
  #--output distortopia/genomes/f1_simulations/rep1_long_read_alignment.vcf.tsv \
  #--rep 1
