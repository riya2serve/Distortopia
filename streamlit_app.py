import streamlit as st
import subprocess
import os
import concurrent.futures
import pandas as pd
import io

from distortopia.simulate_long_reads import simulate_long_reads
from distortopia.align_reads import run_alignment
from distortopia.variant_calling import call_variants
from distortopia.simulate_f1 import generate_f1_from_files as generate_f1_hybrid
from distortopia.snp_detection import run_pipeline
from distortopia.f1_variant_calling import call_f1_variants
from distortopia.compare_variants import load_vcf_as_df, generate_snp_marker_table
from distortopia.detect_crossovers import detect_crossovers

st.set_page_config(page_title="Distortopia: Simulate F1 Genome", layout="centered")
st.title("🌱 Distortopia: Simulate F1 Hybrid and Map Recombination")

# --- Simulate long reads ---
st.subheader("1️⃣ Simulate Long Reads")

col1, col2 = st.columns(2)

with col1:
    if st.button("📡 Simulate A_lyrata Reads"):
        with st.spinner("Simulating PacBio-style reads from A_lyrata..."):
            out1 = simulate_long_reads("A_lyrata.fna", "sim_lyrata.fq.gz", coverage="60x")
            st.success("✅ A_lyrata reads generated.")
            with open(out1, "rb") as f:
                st.download_button("⬇️ Download sim_lyrata.fq.gz", f, file_name="sim_lyrata.fq.gz")

with col2:
    if st.button("📡 Simulate A_thaliana Reads"):
        with st.spinner("Simulating PacBio-style reads from A_thaliana..."):
            out2 = simulate_long_reads("A_thaliana.fna", "sim_thaliana.fq.gz", coverage="60x")
            st.success("✅ A_thaliana reads generated.")
            with open(out2, "rb") as f:
                st.download_button("⬇️ Download sim_thaliana.fq.gz", f, file_name="sim_thaliana.fq.gz")

# --- Check for required input files ---
required_refs = {"A_lyrata.fna", "A_thaliana.fna"}
required_fqs = {"sim_lyrata.fq.gz", "sim_thaliana.fq.gz"}
available_files = set(os.listdir())

missing_refs = required_refs - available_files
missing_fqs = required_fqs - available_files

if missing_refs:
    st.error(f"❌ Missing reference FASTA files: {missing_refs}")
    st.stop()
if missing_fqs:
    st.error(f"❌ Missing FASTQ files: {missing_fqs}")
    st.stop()

# --- Display input files found ---
st.subheader("📂 Input Files Detected")
for genome in ["lyrata", "thaliana"]:
    st.success(f"🧬 {genome.capitalize()} Reference: A_{genome}.fna")
    st.success(f"📄 Simulated Reads: sim_{genome}.fq.gz")

# --- Alignment + Variant Calling ---
st.subheader("2️⃣ Align Reads and Call SNPs")

def full_pipeline(ref, fq, prefix):
    run_alignment(fq, ref, prefix)
    vcf_out = f"{prefix}.vcf"
    call_variants(ref, f"{prefix}.sort.bam", vcf_out)
    return vcf_out

inputs = [
    ("A_lyrata.fna", "sim_lyrata.fq.gz", "sim_lyrata"),
    ("A_thaliana.fna", "sim_thaliana.fq.gz", "sim_thaliana")
]

if st.button("🚀 Run Alignment + Variant Calling for Both"):
    st.info("Running full pipelines in parallel...")
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(full_pipeline, *args) for args in inputs]
        for future in concurrent.futures.as_completed(futures):
            vcf = future.result()
            if os.path.exists(vcf):
                st.success(f"✅ Generated {vcf}")
                with open(vcf) as f:
                    st.download_button(f"⬇️ Download {vcf}", f.read(), file_name=vcf)
            else:
                st.error(f"❌ Failed to generate {vcf}")

# --- Generate F1 hybrid ---
st.subheader("3️⃣ Simulate F1 Hybrid Genome")

if all(os.path.exists(f) for f in ["A_thaliana.fna", "sim_thaliana.vcf", "sim_lyrata.vcf"]):
    if st.button("🔬 Generate F1 Hybrid (IUPAC encoding)"):
        generate_f1_hybrid("A_thaliana.fna", "sim_thaliana.vcf", "sim_lyrata.vcf", "F1_hybrid2.fna")
        st.success("✅ F1 hybrid generated: F1_hybrid2.fna")
        
        # Prepare fresh BytesIO buffer
        with open("F1_hybrid.fna", "rb") as f:
            fasta_bytes = io.BytesIO(f.read())
        st.download_button("⬇️ Download F1_hybrid.fna", fasta_bytes, file_name="F1_hybrid.fna")

# --- Simulate F1 Long Reads ---
st.subheader("4️⃣ Simulate F1 Hybrid Long Reads")

if os.path.exists("F1_hybrid.fna"):
    if st.button("📡 Simulate Long Reads from F1 Hybrid"):
        with st.spinner("Running Badread and streaming output..."):
            simulate_long_reads("F1_hybrid.fna", "F1_hybrid.fq")
            result = subprocess.run(
                [
                    "python", "distortopia/simulate_hybrid_reads.py",
                    "F1_hybrid.fna", "F1_hybrid"
                ],
                capture_output=True, text=True
            )

        if result.returncode == 0:
            st.success("✅ F1 long reads simulated and gzipped.")
            with open("F1_hybrid.fq.gz", "rb") as f:
                st.download_button("⬇️ Download F1_hybrid.fq.gz", f.read(), file_name="F1_hybrid.fq.gz")
        else:
            st.error("❌ Simulation failed.")
            st.text(result.stderr)
else:
    st.warning("Please generate F1_hybrid.fna before simulating long reads.")

# --- Align F1 Hybrid Reads to Both Parents ---
st.subheader("5️⃣ Align F1 Hybrid Reads to Both Parents")

required_files = ["F1_hybrid.fq.gz", "A_thaliana.fna", "A_lyrata.fna"]
if all(os.path.exists(f) for f in required_files):

    with st.expander("🧬 Run Minimap2 + Samtools for F1 → Parents"):
        if st.button("🔗 Align F1 Reads to A. thaliana and A. lyrata"):
            log_box = st.empty()
            try:
                from distortopia.align_hybrid_reads import align_reads

                align_reads("F1_hybrid.fq.gz", "A_thaliana.fna", "F1_to_thaliana", log_box)
                align_reads("F1_hybrid.fq.gz", "A_lyrata.fna", "F1_to_lyrata", log_box)

                st.success("✅ F1 reads successfully aligned to both parental genomes.")
            except Exception as e:
                st.error("❌ Alignment failed. See logs for details.")
else:
    st.warning("⚠️ Missing one or more required files (F1_hybrid.fq.gz, A_thaliana.fna, A_lyrata.fna)")

# --- Call Variants from F1 Alignments ---
st.subheader("6️⃣ Variant Calling on F1 Alignments")

required_bams = ["F1_to_thaliana.sort.bam", "F1_to_lyrata.sort.bam"]
required_refs = ["A_thaliana.fna", "A_lyrata.fna"]

if all(os.path.exists(f) for f in required_bams + required_refs):
    with st.expander("🧬 Run bcftools on F1 Alignments"):
        if st.button("🧪 Call Variants from BAM files"):
            log_box = st.empty()
            try:
                from distortopia.variant_calling import call_variants
                
                call_f1_variants("F1_to_thaliana.sort.bam", "A_thaliana.fna", "F1_to_thaliana.vcf", log_box)
                call_f1_variants("F1_to_lyrata.sort.bam", "A_lyrata.fna", "F1_to_lyrata.vcf", log_box)
                
                st.success("✅ Variant calling completed for both alignments.")
            except Exception as e:
                st.error(f"❌ Variant calling failed: {e}")
else:
    st.warning("⚠️ Missing files for variant calling (check for .bam and .fna files).")

# --- Compare Variants in F1 to both Parents ---
st.subheader("7️⃣ Step 7: Compare F1 Variants to Parental SNPs")

if all(os.path.exists(f) for f in ["sim_thaliana.vcf", "sim_lyrata.vcf", "F1_to_thaliana.vcf", "F1_to_lyrata.vcf"]):

    # Load VCFs into DataFrames
    f1_thal_df = load_vcf_as_df("F1_to_thaliana.vcf", ref_name="Thaliana", alt_name="F1")
    sim_thal_df = load_vcf_as_df("sim_thaliana.vcf", ref_name="Thaliana", alt_name="Sim")
    f1_lyr_df  = load_vcf_as_df("F1_to_lyrata.vcf", ref_name="Lyrata", alt_name="F1")
    sim_lyr_df = load_vcf_as_df("sim_lyrata.vcf", ref_name="Lyrata", alt_name="Sim")

    # 🔍 Compare F1 to A. thaliana
    st.markdown("#### 🔍 Compare F1 to A. thaliana")
    merged_thal = pd.merge(f1_thal_df, sim_thal_df, on=["CHROM", "POS"], how="outer", indicator=True, suffixes=("_F1", "_Thaliana"))
    shared_thal = merged_thal[merged_thal["_merge"] == "both"]
    f1_unique_thal = merged_thal[merged_thal["_merge"] == "left_only"]

    st.write(f"Shared SNPs: {len(shared_thal)}")
    st.write(f"F1-unique SNPs: {len(f1_unique_thal)}")

    # 🔍 Compare F1 to A. lyrata
    st.markdown("#### 🔍 Compare F1 to A. lyrata")
    merged_lyr = pd.merge(f1_lyr_df, sim_lyr_df, on=["CHROM", "POS"], how="outer", indicator=True, suffixes=("_F1", "_Lyrata"))
    shared_lyr = merged_lyr[merged_lyr["_merge"] == "both"]
    f1_unique_lyr = merged_lyr[merged_lyr["_merge"] == "left_only"]

    st.write(f"Shared SNPs: {len(shared_lyr)}")
    st.write(f"F1-unique SNPs: {len(f1_unique_lyr)}")

    # 🧪 Optional: Write SNP marker table to disk for downstream use
    if not os.path.exists("compare_variants_output.tsv"):
        try:
            marker_df = pd.merge(sim_thal_df, sim_lyr_df, on=["CHROM", "POS"], how="inner")
            marker_df = marker_df[["CHROM", "POS", "Thaliana", "Lyrata"]]
            marker_df.to_csv("compare_variants_output.tsv", sep="\t", index=False)
            st.success("✅ SNP marker table saved as compare_variants_output.tsv")
        except Exception as e:
            st.error(f"❌ Failed to write SNP marker table: {e}")
else:
    st.warning("❌ Required VCF files not found. Please complete earlier steps first.")

# --- Detect Crossovers ---
st.subheader("8️⃣ Step 8: Detect Crossovers from F1 Reads")

if os.path.exists("compare_variants_output.tsv"):
    
    tab1, tab2 = st.tabs(["🧬 F1 → A. thaliana", "🧬 F1 → A. lyrata"])

    with tab1:
        if os.path.exists("F1_to_thaliana.sort.bam"):
            if st.button("Run Crossover Detection (Thaliana)"):
                with st.spinner("Running crossover detection on F1_to_thaliana.sort.bam..."):
                    detect_crossovers(
                        bam_path="F1_to_thaliana.sort.bam",
                        marker_table_path="compare_variants_output.tsv",
                        output_path="crossovers_thaliana.tsv"
                    )
                st.success("✅ Thaliana crossover detection complete. Results saved to `crossovers_thaliana.tsv`.")

            if os.path.exists("crossovers_thaliana.tsv"):
                co_df_thal = pd.read_csv("crossovers_thaliana.tsv", sep="\t")
                st.write(f"🔍 Detected {co_df_thal['num_crossovers'].sum()} crossovers in {len(co_df_thal)} reads.")
                if st.checkbox("Show Thaliana crossover calls"):
                    st.dataframe(co_df_thal[["qname", "chrom", "start", "end", "pattern", "num_crossovers"]].head(100))
        else:
            st.warning("❌ File `F1_to_thaliana.sort.bam` not found.")

    with tab2:
        if os.path.exists("F1_to_lyrata.sort.bam"):
            if st.button("Run Crossover Detection (Lyrata)"):
                with st.spinner("Running crossover detection on F1_to_lyrata.sort.bam..."):
                    detect_crossovers(
                        bam_path="F1_to_lyrata.sort.bam",
                        marker_table_path="compare_variants_output.tsv",
                        output_path="crossovers_lyrata.tsv"
                    )
                st.success("✅ Lyrata crossover detection complete. Results saved to `crossovers_lyrata.tsv`.")

            if os.path.exists("crossovers_lyrata.tsv"):
                co_df_lyr = pd.read_csv("crossovers_lyrata.tsv", sep="\t")
                st.write(f"🔍 Detected {co_df_lyr['num_crossovers'].sum()} crossovers in {len(co_df_lyr)} reads.")
                if st.checkbox("Show Lyrata crossover calls"):
                    st.dataframe(co_df_lyr[["qname", "chrom", "start", "end", "pattern", "num_crossovers"]].head(100))
        else:
            st.warning("❌ File `F1_to_lyrata.sort.bam` not found.")

else:
    st.warning("❌ SNP marker table `compare_variants_output.tsv` not found. Run previous steps first.")

# --- Visualize ALL VCF tables ---
st.subheader("🧬 SNP Tables from Variant Calling")

vcf_files = ["sim_lyrata.vcf", 
             "sim_thaliana.vcf",
             "F1_to_lyrata.vcf",
             "F1_to_thaliana.vcf"
]

for vcf_file in vcf_files:
    if os.path.exists(vcf_file):  
        st.markdown(f"**{vcf_file}**")

        # Let user choose how many rows to preview
        num_rows = st.slider(f"Rows to preview from {vcf_file}", 50, 100, 500, step=100)

        with open(vcf_file) as f:
            lines = [l for l in f if not l.startswith("##")]

        df = pd.read_csv(io.StringIO("".join(lines)), sep="\t", nrows=num_rows)
        st.dataframe(df)
        st.caption(f"Showing first {num_rows} rows")

