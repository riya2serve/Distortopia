import os
import argparse
import zipfile
import shutil
from ncbi.datasets import GenomeApi
from Bio import Entrez

# ==================
# CONFIGURATION
# ==================

Entrez.email = "rr3491@columbia.edu"  # Required by NCBI
output_folder = "user-data"
os.makedirs(output_folder, exist_ok=True)

# ==================
# FUNCTIONS
# ==================

def get_species():
    """
    Prompts user for species names via terminal.
    """
    species_dict = {}
    print("Enter plant species name to search on NCBI (e.g. 'Amaranthus retroflexus'). Type 'done' to finish.")

    while True:
        species = input("Species name (or 'done'): ").strip()
        if species.lower() == "done":
            break
        species_dict[species] = None
    return species_dict
   
def fetch_genomes(species_name, output_folder="user-data", force_download=False):
    """
    Uses the NCBI Datasets API to fetch genome FASTA and GFF files as a ZIP.
    """
    from ncbi.datasets import V1alpha1AssemblyCatalogApi, GenomeApi
    os.makedirs(output_folder, exist_ok=True)
    safe_name = species_name.replace(" ", "_")
    output_zip = os.path.join(output_folder, f"{safe_name}.zip")

    if os.path.exists(output_zip) and not force_download:
        print(f"{output_zip} already exists. Use --force to re-download.")
        return

    try:
        print(f"Searching for best assembly for {species_name}...")
        catalog_api = V1alpha1AssemblyCatalogApi()
        search_result = catalog_api.assembly_catalog_search(taxon=species_name, page_size=1)
        assemblies = search_result.assemblies

        if not assemblies:
            print(f"No genome assemblies found for {species_name}")
            return

        accession = assemblies[0].assembly.accession
        print(f"Found accession {accession} for {species_name}")

        api = GenomeApi()
        print(f"Downloading genome for {species_name} using accession {accession}...")

        response = api.download_assembly_package(
            accessions=[accession],
            include_annotation_type=["GENOME_FASTA", "GENOME_GFF"],
            filename=output_zip,
            format="zip"
        )

        with open(output_zip, "wb") as f:
            shutil.copyfileobj(response, f)

        print(f"Successfully downloaded genome for {species_name} → {output_zip}")

    except Exception as e:
        print(f"Failed to download genome for {species_name}: {e}")

def extract_genomes(zip_path, extract_to="user-data"):
    """
    Optional: Extracts downloaded zip files into the output directory.
    """
    try:
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(extract_to)
        print(f"Extracted contents of {zip_path}")
    except Exception as e:
        print(f"Failed to extract {zip_path}: {e}")

def parse_args():
    """
    Handles command-line arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true", help="Force re-download even if files exist")
    parser.add_argument("--extract", action="store_true", help="Extract downloaded ZIP files")
    return parser.parse_args()

# ==================
# MAIN PROGRAM
# ==================

if __name__ == "__main__":
    args = parse_args()
    species_dict = get_species()

    if not species_dict:
        print("No species provided. Exiting.")
    else:
        for species in species_dict.keys():
            fetch_genomes(species, output_folder=output_folder, force_download=args.force)
            if args.extract:
                zip_path = os.path.join(output_folder, f"{species.replace(' ', '_')}.zip")
                extract_genomes(zip_path, extract_to=output_folder)

        print("All downloads finished.")






