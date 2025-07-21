# distortopia/simulate_hybrid_reads.py

# distortopia/simulate_hybrid_reads.py

import subprocess
import argparse
import gzip
import shutil
import os
from datetime import datetime

def log(message, log_path="simulate_hybrid_reads.log"):
    timestamp = datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")
    with open(log_path, "a") as f:
        f.write(f"{timestamp} {message}\n")

def simulate_long_reads(reference_fasta, output_fq, coverage="60x",
                        read_length_mean=8000, read_length_sd=3000,
                        error_model="pacbio2021", gzip_output=True):

    cmd = [
        "badread", "simulate",
        "--reference", reference_fasta,
        "--quantity", coverage,
        "--length", f"{read_length_mean},{read_length_sd}",
        "--error_model", error_model
    ]

    # Streamlit placeholder for real-time output
    log_box = st.empty()
    log_box.text(f"🚀 Running Badread: {' '.join(cmd)}")

    # Start process
    try:
        with open(output_fq, "w") as f_out:
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

            for line in proc.stdout:
                f_out.write(line)
                log_box.text(line.strip())  # update live in Streamlit

            proc.wait()
            if proc.returncode != 0:
                log_box.text("❌ Badread simulation failed.")
                raise RuntimeError("Badread simulation failed.")
        
        log_box.text(f"✅ Badread wrote: {output_fq}")

        # Compress if needed
        if gzip_output:
            out_fq_gz = output_fq + ".gz"
            with open(output_fq, "rb") as f_in, gzip.open(out_fq_gz, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.remove(output_fq)
            log_box.text(f"✅ Output gzipped: {out_fq_gz}")
            return out_fq_gz

        return output_fq

    except Exception as e:
        log_box.text(f"❌ Error: {e}")
        raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Simulate F1 hybrid long reads with Badread")
    parser.add_argument("reference", help="Path to F1_hybrid.fna")
    parser.add_argument("output_prefix", help="Prefix for output (e.g., F1_hybrid)")
    args = parser.parse_args()

    simulate_long_reads(args.reference, args.output_prefix)
    print(f"✅ F1 long reads written to {args.output_prefix}.fq.gz")

