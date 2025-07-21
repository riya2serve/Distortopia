import subprocess

def call_variants(reference_fasta, sorted_bam, output_vcf):
    """
    Call high-confidence SNPs using bcftools with multithreading and filtering.
    """
    print(f"[INFO] Calling variants from {sorted_bam} → {output_vcf}")

    with open(output_vcf, "w") as vcf_out:
        mpileup = subprocess.Popen(
            ["bcftools", "mpileup", "-Ou", "-f", reference_fasta,
             "--threads", "8", "--max-depth", "250", sorted_bam],
            stdout=subprocess.PIPE
        )
        call = subprocess.Popen(
            ["bcftools", "call", "-mv", "-Ou", "--threads", "8"],
            stdin=mpileup.stdout,
            stdout=subprocess.PIPE
        )
        view = subprocess.run(
            ["bcftools", "view", "-Ov", "-v", "snps", "-m2", "-M2",
             "-i", "QUAL>20 && DP>10"],
            stdin=call.stdout,
            stdout=vcf_out,
            check=True
        )

    print(f"[DONE] Wrote filtered VCF: {output_vcf}")

