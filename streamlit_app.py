import streamlit as st
import subprocess
import os
import concurrent.futures
import pandas as pd
import io
import matplotlib.pyplot as plt

st.set_page_config(page_title="Distortopia Variant Simulation", layout="centered")
st.title("Distortopia: Simulate F1 Hybrid Genome")

# --- Simulate long reads ---
from distortopia.simulate_long_reads import simulate_long_reads

st.subheader("Simulate Long Reads from Reference Genomes")

col1, col2 = st.columns(2)

with col1:
    if st.button("Simulate A_lyrata Reads"):
        with st.spinner("Simulating A_lyrata long reads..."):
            out1 = simulate_long_reads("A_lyrata.fna", "sim_lyrata.fq.gz", coverage="60x")
            st.success("A_lyrata reads generated.")
            with open(out1, "rb") as f:
                st.download_button("Download A_lyrata.fq.gz", f, file_name="sim_lyrata.fq.gz")

with col2:
    if st.button("Simulate A_thaliana Reads"):
        with st.spinner("Simulating A_thaliana long reads..."):
            out2 = simulate_long_reads("A_thaliana.fna", "sim_thaliana.fq.gz", coverage="60x")
            st.success("A_thaliana reads generated.")
            with open(out2, "rb") as f:
                st.download_button("Download A_thaliana.fq.gz", f, file_name="sim_thaliana.fq.gz")

# --- Detect input files ---
fnas = sorted([f for f in os.listdir() if f.endswith(".fna")])
fq_files = os.listdir()
fqs = sorted([
    f for f in fq_files
    if f.endswith(".fq") or f.endswith("fq.gz")
])

required_refs = {"A_lyrata.fna", "A_thaliana.fna"}
required_fqs = {"sim_lyrata.fq.gz", "sim_thaliana.fq.gz"}

missing_fqs = required_fqs - set(fqs)

if not required_refs.issubset(set(fnas)):
    st.error(f"Missing required reference FASTA files: {required_refs - set(fnas)}")
    st.stop()

if missing_fqs:
    st.error(f"Missing required FASTQ files: {missing_fqs}")
    st.stop()

# --- Show detected files ---
inputs = [
    ("A_lyrata.fna", "sim_lyrata.fq.gz", "sim_lyrata.vcf"),
    ("A_thaliana.fna", "sim_thaliana.fq.gz", "sim_thaliana.vcf")
]

for ref, fq, _ in inputs:
    st.success(f"Found reference genome: {ref}")
    st.success(f"Found FASTQ file: {fq}")

# --- Variant calling function ---
def run_pipeline(ref, fq, vcf_out):
    sam = fq.replace(".fq", ".sam")
    bam = fq.replace(".fq", ".bam")
    sorted_bam = fq.replace(".fq", ".sort.bam")

    subprocess.run(["minimap2", "-a", ref, fq], stdout=open(sam, "w"), check=True)
    subprocess.run(["samtools", "view", "-bS", sam], stdout=open(bam, "wb"), check=True)
    subprocess.run(["samtools", "sort", "-o", sorted_bam, bam], check=True)
    subprocess.run(["samtools", "index", sorted_bam], check=True)

    with open(vcf_out, "w") as vcf:
        mpileup = subprocess.Popen(["bcftools", "mpileup", "-f", ref, sorted_bam], stdout=subprocess.PIPE)
        subprocess.run(["bcftools", "call", "-mv", "-Ov"], stdin=mpileup.stdout, stdout=vcf, check=True)

    return vcf_out

# --- Run pipelines ---
if st.button("Run both pipelines"):
    st.info("Running alignments and variant calling in parallel...")
    vcf_files = []
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(run_pipeline, *args) for args in inputs]
        for future in concurrent.futures.as_completed(futures):
            vcf_file = future.result()
            vcf_files.append(vcf_file)
            if os.path.exists(vcf_file):
                st.success(f" Completed: {vcf_file}")
                with open(vcf_file) as f:
                    st.download_button(f"Download {vcf_file}", f.read(), file_name=vcf_file)
            else:
                st.error(f" Failed to generate {vcf_file}")

# --- Generate F1 Hybrid (dual-parent) ---
from generate_f1_hybrid import generate_f1_hybrid

if all(os.path.exists(f) for f in ["A_thaliana.fna", "sim_thaliana.vcf", "sim_lyrata.vcf"]):
    if st.button("Generate F1 Hybrid Genome (dual-parent with IUPAC)"):
        generate_f1_hybrid("A_thaliana.fna", "sim_thaliana.vcf", "sim_lyrata.vcf", "F1_hybrid.fna")
        st.success("F1 hybrid genome simulated as F1_hybrid.fna")

        with open("F1_hybrid.fna") as f:
            st.download_button("Download F1_hybrid.fna", f.read(), file_name="F1_hybrid.fna")

# --- Display VCF Tables ---
for vcf in ["sim_lyrata.vcf", "sim_thaliana.vcf"]:
    if os.path.exists(vcf):
        st.subheader(f"Variant Table: {vcf}")
        with open(vcf) as f:
            lines = [line for line in f if not line.startswith("##")]
        df = pd.read_csv(io.StringIO("".join(lines)), sep='\t')
        st.dataframe(df)



