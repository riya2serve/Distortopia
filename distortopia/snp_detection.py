import subprocess
import os

def run_pipeline(ref, fq, vcf_out):
    sam = fq.replace(".fq", ".sam")
    bam = fq.replace(".fq", ".bam")
    sorted_bam = fq.replace(".fq", ".sort.bam")

    # Ensure FASTA index exists
    if not os.path.exists(f"{ref}.fai"):
        subprocess.run(["samtools", "faidx", ref], check=True)

    # Align and convert to sorted BAM
    subprocess.run(["minimap2", "-ax", "map-pb", ref, fq], stdout=open(sam, "w"), check=True)
    subprocess.run(["samtools", "view", "-bS", sam], stdout=open(bam, "wb"), check=True)
    subprocess.run(["samtools", "sort", "-o", sorted_bam, bam], check=True)
    subprocess.run(["samtools", "index", sorted_bam], check=True)

    # Variant calling with relaxed thresholds
    with open(vcf_out, "w") as vcf:
        mpileup = subprocess.Popen([
            "bcftools", "mpileup",
            "-f", ref,
            "--max-depth", "10000",
            "-a", "AD,DP",
            sorted_bam
        ], stdout=subprocess.PIPE)

        subprocess.run([
            "bcftools", "call",
            "-mv",    # call both SNPs and indels
            "-Ov"     # output VCF in text format
        ], stdin=mpileup.stdout, stdout=vcf, check=True)

    return vcf_out

