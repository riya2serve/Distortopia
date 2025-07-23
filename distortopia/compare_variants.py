import pandas as pd
import io

def load_vcf_as_df(vcf_path, ref_name="REF", alt_name="ALT"):
    with open(vcf_path) as f:
        lines = [l for l in f if not l.startswith("##")]
    df = pd.read_csv(io.StringIO("".join(lines)), sep="\t")
    df = df[["#CHROM", "POS", "REF", "ALT"]]
    df.rename(columns={
        "#CHROM": "CHROM",
        "REF": ref_name,
        "ALT": alt_name
    }, inplace=True)
    return df

def generate_snp_marker_table(thaliana_vcf, lyrata_vcf, outname="snp_marker_table.tsv"):
    df_thal = load_vcf_as_df(thaliana_vcf, "Thaliana_REF", "Thaliana_ALT")
    df_lyra = load_vcf_as_df(lyrata_vcf, "Lyrata_REF", "Lyrata_ALT")

    print("🔍 Thal SNPs:", len(df_thal))
    print("🔍 Lyra SNPs:", len(df_lyra))

    merged = pd.merge(df_thal, df_lyra, on=["CHROM", "POS"], how="inner")
    print("✅ Shared SNPs:", len(merged))

    marker_table = merged[["CHROM", "POS", "Thaliana_REF", "Lyrata_ALT"]]
    marker_table.rename(columns={
        "Thaliana_REF": "Thaliana",
        "Lyrata_ALT": "Lyrata"
    }, inplace=True)

    marker_table.to_csv(outname, sep="\t", index=False)
    print(f"📁 Saved SNP marker table: {outname}")

    return marker_table

if __name__ == "__main__":
    thaliana_vcf = "sim_thaliana.vcf"
    lyrata_vcf = "sim_lyrata.vcf"
    generate_snp_marker_table(thaliana_vcf, lyrata_vcf)

