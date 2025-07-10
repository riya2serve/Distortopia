import streamlit as st
import subprocess
import os
import concurrent.futures
import pandas as pd
import io
import matplotlib.pyplot as plt
import recombination_rate

st.set_page_config(page_title="COmapper Variant Calling", layout="centered")
st.title("COmapper: Local Variant Calling Test")

# --- Detect input files ---
fnas = sorted([f for f in os.listdir() if f.endswith(".fna")])
fqs = sorted([f for f in os.listdir() if f.endswith(".fq")])

# --- Require specific files ---
required_refs = {"A_lyrata.fna", "A_thaliana.fna"}
required_fqs = {"sim_lyrata.fq", "sim_thaliana.fq"}

if not required_refs.issubset(fnas) or not required_fqs.issubset(fqs):
    st.error("Missing one or more required files: A_lyrata.fna, A_thaliana.fna, sim_lyrata.fq, sim_thaliana.fq")
    st.stop()

# --- Input combinations ---
inputs = [
    ("A_lyrata.fna", "sim_lyrata.fq", "sim_lyrata.vcf"),
    ("A_thaliana.fna", "sim_thaliana.fq", "sim_thaliana.vcf")
]

# --- Display detected files ---
for ref, fq, _ in inputs:
    st.success(f"Found reference genome: {ref}")
    st.success(f"Found FASTQ file: {fq}")

# --- Variant calling pipeline function ---
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

# --- Run pipelines when button is clicked ---
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

# --- Optional: Display most recent VCF ---
for vcf in ["sim_lyrata.vcf", "sim_thaliana.vcf"]:
    if os.path.exists(vcf):
        st.subheader(f"Variant Table: {vcf}")
        with open(vcf) as f:
            lines = [line for line in f if not line.startswith("##")]
        df = pd.read_csv(io.StringIO("".join(lines)), sep='\t')
        st.dataframe(df)

# --- Load all 3 VCFs ---
def load_vcf_positions(vcf_file):
    with open(vcf_file) as f:
        lines = [l for l in f if not l.startswith("##")]
    if not lines:
        return set()
    df = pd.read_csv(io.StringIO("".join(lines)), sep="\t")
    return set(zip(df["#CHROM"], df["POS"]))

if st.checkbox("Compare variant positions across F1 and parents"):

    f1_pos = load_vcf_positions("F1.vcf")
    thaliana_pos = load_vcf_positions("sim_thaliana.vcf")
    lyrata_pos = load_vcf_positions("sim_lyrata.vcf")

    all_pos = {
        "F1": f1_pos,
        "Thaliana": thaliana_pos,
        "Lyrata": lyrata_pos
    }

    shared_all = f1_pos & thaliana_pos & lyrata_pos
    f1_unique = f1_pos - (thaliana_pos | lyrata_pos)
    thaliana_unique = thaliana_pos - (f1_pos | lyrata_pos)
    lyrata_unique = lyrata_pos - (f1_pos | thaliana_pos)

    st.write(f" Shared in all 3: {len(shared_all)}")
    st.write(f" F1 only: {len(f1_unique)}")
    st.write(f" Thaliana only: {len(thaliana_unique)}")
    st.write(f" Lyrata only: {len(lyrata_unique)}")

    # --- Visualization ---
    labels = ["Shared All", "F1 only", "Thaliana only", "Lyrata only"]
    counts = [len(shared_all), len(f1_unique), len(thaliana_unique), len(lyrata_unique)]

    fig, ax = plt.subplots()
    ax.bar(labels, counts)
    ax.set_ylabel("SNP count")
    ax.set_title("Shared vs. Unique Variant Positions")
    st.pyplot(fig)

if st.checkbox("Estimate recombination rate from VCFs"):
    switches, informative_sites = recombination_rate.estimate_recombination_rate_v2(
        "sim_thaliana.vcf", "sim_lyrata.vcf", "F1.vcf"
    )
    st.subheader("Recombination Summary")
    st.write(f"Number of informative sites: {informative_sites}")
    st.write(f"Number of allele switches (crossovers): {switches}")
    if informative_sites > 1:
        rate = switches / informative_sites
        st.write(f"Estimated recombination rate: **{rate:.4f} crossovers per site**")


