import subprocess
import os

def run_alignment_pipeline(name, ref_path, fq_gz_path, workdir=".", corenum=4):
    """
    Full COmapper-style alignment and filtering pipeline.

    Args:
        name (str): Dataset name (e.g., 'sim_lyrata')
        ref_path (str): Path to reference genome index (.mmi)
        fq_gz_path (str): Path to .fq.gz input file
        workdir (str): Output directory
        corenum (int): Number of threads to use
    """

    os.makedirs(workdir, exist_ok=True)
    os.makedirs(f"{workdir}/nanoplot", exist_ok=True)
    os.makedirs(f"{workdir}/Q14_output", exist_ok=True)

    fq_path = os.path.join(workdir, f"{name}.fq")
    fq_copy = os.path.join(workdir, f"{name}.fq.gz")
    sam_path = os.path.join(workdir, f"{name}.sam")
    bam_path = os.path.join(workdir, f"{name}.bam")
    sorted_bam = os.path.join(workdir, f"{name}.sort.bam")
    sbb_bam = os.path.join(workdir, f"{name}.sort.sbb.bam")
    sbb_sam = os.path.join(workdir, f"{name}.sort.sbb.sam")

    # Copy and unzip .fq.gz
    subprocess.run(["cp", fq_gz_path, fq_copy], check=True)
    subprocess.run(["gunzip", fq_copy], check=True)
    print("Unzipped FASTQ")

    # Alignment
    subprocess.run([
        "minimap2", "-t", str(corenum), "-ax", "map-ont", "--cs",
        ref_path, fq_path
    ], stdout=open(sam_path, "w"), check=True)
    print("Alignment done")

    # Convert to BAM
    subprocess.run(["samtools", "view", "-@", str(corenum), "-Sb"]()

