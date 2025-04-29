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
    Returns a dictionary of species names. 
    """
    species_dict = {}
    print("Enter plant species name to search on NCBI. Type 'done' to finish.")

    while True:
        species = input("Species name (or 'done'): ").strip()
        if species.lower() == "done":
            break
        species_dict[species] = None  # For storing accessions later
    return species_dict

def fetch_genomes(species_name, output_folder="user-data", force_download=False, extract=False):
    """
    Downloads a dehydrated .zip file containing metadata and fetch.txt.
    Then rehydrates it to pull actual genome and GFF files.
    Optionally extracts .fna.gz and .gff.gz from the final output.
    """
    # Sanitize species name for folder/filename use
    safe_name = species_name.replace(" ", "_")
    zip_path = os.path.join(output_folder, f"{safe_name}.zip")
    unzip_dir = os.path.join(output_folder, safe_name)

    # If zip already exists and --force not specified, skip download
    if os.path.exists(zip_path) and not force_download:
        print(f"{zip_path} already exists. Use --force to re-download.")
    else:
        print(f"Downloading dehydrated genome for: {species_name}")

        try:
            # Step 1: Download a dehydrated .zip (contains metadata & fetch.txt)
            subprocess.run([
                "datasets", "download", "genome", "taxon", species_name,
                "--filename", zip_path,
                "--dehydrated",
                "--include", "genome,gff3"
            ], check=True)
            print(f"Downloaded dehydrated archive to: {zip_path}")
        except subprocess.CalledProcessError as e:
            print(f"Download failed for {species_name}: {e}")
            return

    try:
        # Step 2: Unzip the dehydrated archive to a working folder
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(unzip_dir)
            print(f"Unzipped to: {unzip_dir}")
    except Exception as e:
        print(f"Failed to unzip: {e}")
        return

    try:
        # Step 3: Rehydrate → uses fetch.txt to download actual genome + annotation files
        subprocess.run(["datasets", "rehydrate", "--directory", unzip_dir], check=True)
        print(f"Rehydrated genome into: {unzip_dir}")
    except subprocess.CalledProcessError as e:
        print(f"Rehydration failed: {e}")
        return

    # Step 4: (Optional) Extract .fna.gz and .gff.gz files from rehydrated structure
    if extract:
        extract_fasta_gff(unzip_dir)

def extract_fasta_gff(base_dir):
    """
    Find and extract .fna.gz and .gff.gz files from rehydrated dataset directory.
    """
    try:
        with zipfile.ZipFile(os.path.join(base_dir, "ncbi_dataset.zip"), 'r') as zip_ref:
            files_to_extract = [f for f in zip_ref.namelist() if f.endswith(('.fna.gz', '.gff.gz'))]
            if files_to_extract:
                print(f"Extracting files: {', '.join(files_to_extract)}")
                zip_ref.extractall(base_dir, members=files_to_extract)
                print(f"Extracted .fna.gz and .gff.gz to: {base_dir}")
            else:
                print("No FASTA or GFF files found in rehydrated ZIP.")
    except Exception as e:
        print(f"Extraction failed: {e}")

def parse_args():
    """
    Parse command line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true", help="Force re-download even if zip exists")
    parser.add_argument("--extract", action="store_true", help="Extract .fna.gz and .gff.gz after rehydration")
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
        for species in species_dict:
            fetch_genomes(
                species,
                output_folder=output_folder,
                force_download=args.force,
                extract=args.extract
            )

        print("All downloads and processing complete.")




