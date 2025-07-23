# distortopia/f1_variant_calling.py

import subprocess
import os

def log(msg, log_box=None):
    """Display message to console and optionally to Streamlit."""
    print(msg)
    if log_box:
        log_box.info(msg)

def call_f1_variants(bam_file, reference_fasta, output_vcf, log_box=None, threads=4):
    """Call variants from F1 alignments using bcftools mpileup and call."""
    log(f"📥 Starting variant calling for {bam_file} vs {reference_fasta}", log_box)

    try:
        # Construct bcftools mpileup command
        mpileup_cmd = [
            "bcftools", "mpileup",
            "--threads", str(threads),
            "-Ou",  # Output uncompressed BCF
            "-f", reference_fasta,
            bam_file
        ]

        # Construct bcftools call command
        call_cmd = [
            "bcftools", "call",
            "--threads", str(threads),
            "-mv",  # Output only variant sites
            "-Ov",  # Output uncompressed VCF
            "-o", output_vcf
        ]

        log("🔬 Running bcftools mpileup + call...", log_box)
        log(f"  🧬 mpileup: {' '.join(mpileup_cmd)}", log_box)
        log(f"  🧪 call:    {' '.join(call_cmd)}", log_box)

        # Run mpileup and pipe to call
        p1 = subprocess.Popen(mpileup_cmd, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(call_cmd, stdin=p1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p1.stdout.close()
        out, err = p2.communicate()

        if p2.returncode != 0:
            log(f"❌ bcftools call failed:\n{err.decode()}", log_box)
            raise RuntimeError("Variant calling failed.")

        if os.path.exists(output_vcf):
            log(f"✅ Variants successfully written to {output_vcf}", log_box)
        else:
            raise RuntimeError("VCF file was not created.")

    except Exception as e:
        log(f"❌ Error during variant calling: {e}", log_box)
        raise

