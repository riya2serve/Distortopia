import os
import argparse
import subprocess
import zipfile

# ================
# CONFIGURATION
# ================

output_folder = "user-data"
os.makedirs(output_folder, exist_ok=True)  # Create the directory if it doesn't exist

# ================
# MY FUNCTIONS
# ================

def get_species():
    """
    Prompt the user to input species names for searching in NCBI.
    """
    species_dict = {}
    print("Enter plant species name to search on NCBI. Type 'done' to finish.")

    while True:
        species = input("Species name (or 'done'): ").strip()
        if species.lower() == "done":
            break
        species_dict[species] = None  # For storing accessions later
    return species_dict

def fetch_genomes(species_name, output_folder="user-data", force_download=False):
    """
    Fetch the genome assembly for the given species using the NCBI Datasets CLI tool.
    Downloads and optionally extracts the files.
    """
    safe_name = species_name.replace(" ", "_")
    zip_path = os.path.join(output_folder, f"{safe_name}.zip")

    if os.path.exists(zip_path) and not force_download:
        print(f"{zip_path} already exists. Use --force to re-download.")
        return

    print(f"Using NCBI Datasets CLI to fetch genome for: {species_name}")

    try:
        subprocess.run([
            "datasets", "download", "genome", "taxon", species_name,
            "--filename", zip_path,
            "--dehydrated",
            "--include", "genome,gff3"
            ], check=True)

        #Rehydrate after download
        subprocess.run([
            "datasets", "rehydrate",
            "--directory", output_folder
            ], check=True)

        print(f"Downloaded genome for {species_name} → {zip_path}")

        if args.extract:
            extract_zip(zip_path)

    except subprocess.CalledProcessError as e:
        print(f"Failed to download genome for {species_name}: {e}")

def extract_zip(zip_path, extract_to="user-data"):
    """
    Extract only the FASTA (.fna.gz) and GFF (.gff.gz) files from the downloaded zip file.
    """
    try:
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            # List of files to extract (FASTA and GFF files only)
            files_to_extract = [f for f in zip_ref.namelist() if f.endswith(('.fna.gz', '.gff.gz'))]

            if files_to_extract:
                print(f"Extracting files: {', '.join(files_to_extract)}")
                for file in files_to_extract:
                    zip_ref.extract(file, extract_to)
                print(f"Extracted specified files to: {extract_to}")
            else:
                print("No FASTA or GFF files found in the zip archive.")
    except Exception as e:
        print(f"Failed to extract ZIP: {e}")

def parse_args():
    """
    Parse command line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true", help="Force re-download even if files exist")
    parser.add_argument("--extract", action="store_true", help="Extract downloaded zip files")
    return parser.parse_args()

# ==================
# MAIN PROGRAM LOGIC
# ==================

if __name__ == "__main__":
    args = parse_args()
    species_dict = get_species()

    if not species_dict:
        print("No species provided. Exiting.")
    else:
        for species in species_dict.keys():
            fetch_genomes(species, output_folder=output_folder, force_download=args.force)

        print("All downloads completed.")





