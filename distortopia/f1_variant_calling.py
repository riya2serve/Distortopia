# distortopia/f1_variant_calling.py

import subprocess
import os

def log(msg, log_box=None):
    print(msg)
    if log_box:
        log_box.info(msg)

def call_f1_variants(bam_file, reference_fasta, output_vcf, log_box=None):
    log(f"📥 Calling variants with bcftools for {bam_file}", log_box)
    # your subprocess calls here...
    
    try:
        # Create a BCF file using bcftools mpileup
        bcf_file = output_vcf.replace(".vcf", ".bcf")
        mpileup_cmd = [
            "bcftools", "mpileup", "-Ou",
            "-f", reference_fasta,
            bam_file
        ]
        call_cmd = [
            "bcftools", "call", "-mv", "-Ov", "-o", output_vcf
        ]

        log("🔬 Running bcftools mpileup + call...", log_box)
        p1 = subprocess.Popen(mpileup_cmd, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(call_cmd, stdin=p1.stdout)
        p1.stdout.close()
        p2.communicate()

        if os.path.exists(output_vcf):
            log(f"✅ Variants written to {output_vcf}", log_box)
        else:
            raise RuntimeError("VCF file not created.")

    except Exception as e:
        log(f"❌ Variant calling failed: {e}", log_box)
        raise RuntimeError("Variant calling failed.") from e

