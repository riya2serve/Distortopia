import subprocess
import os

def simulate_long_reads(reference_fasta, output_fq, coverage="60x", 
                        read_length_mean=8000, read_length_sd=3000,
                        error_model="pacbio2021", gzip_output=True):
    """
    Runs Badread to simulate long reads from a reference FASTA.
    
    Args:
        reference_fasta (str): Path to reference .fna file
        output_fq (str): Output FASTQ file path (e.g., sim_lyrata.fq)
        coverage (str): e.g., "60x"
        read_length_mean (int): Mean read length
        read_length_sd (int): Standard deviation of read length
        error_model (str): Badread error model (e.g., "pacbio2021")
        gzip_output (bool): Whether to gzip the output
    """
    
    cmd = [
        "badread", "simulate",
        "--reference", reference_fasta,
        "--quantity", coverage,
        "--length", f"{read_length_mean},{read_length_sd}",
        "--error_model", error_model
    ]

    # Open gzip pipe if requested
    if gzip_output:
        out_fq_gz = output_fq + ".gz"
        with open(out_fq_gz, "wb") as f:
            p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            p2 = subprocess.Popen(["gzip"], stdin=p1.stdout, stdout=f)
            p1.stdout.close()
            p2.communicate()
        return out_fq_gz
    else:
        with open(output_fq, "w") as f:
            subprocess.run(cmd, stdout=f, check=True)
        return output_fq
