import os 
import argparse
from ftplib import FTP
import gzip #useful for FASTA and GFF file decompression later on 
from Bio import Entrez #for NCBI access
from Bio import SeqIO #to parse through FASTA files later on

# ================
# CONFIGURATION
# ================

Entrez.email = "rr3491@columbia.edu" #required by NCBI 

# ================
# MY FUNCTIONS
# ================

#get species input from user via terminal:
def get_species():
    species_dict = {}
    print("Enter plant species name to search on NCBI (e.g. 'Amaranthus retrofleuxs'). Type 'done' to finish.")

    while True:
        species = input("Species name (or 'done'): ").strip()
        if species.lower() == "done":
            break
        species_dict[species] = None #for accessions 
    return species_dict 

#use Entrez module to fetch genome assembly accession from NCBI
def get_assembly(species_name):
    """
    Get lastest assembly accessions for each given species name.
    """
    try:
        handle = Entrez.esearch(db = "assembly", term = species_name + "[Organism]", retmax = 1, sort = "relevance")
        record = Entrez.read(handle)
        handle.close()
        if not record["IdList"]:
            print(f"[ERROR] No assemblies found for {species_name}")
            return None

        uid = record["IdList"][0]

        handle = Entrez.esummary(db="assembly", id=uid)
        summary = Entrez.read(handle)
        handle.close()

        doc = summary['DocumentSummarySet']['DocumentSummary'][0]
        acc = doc['AssemblyAccession']
        print(f"[INFO] Found assembly accession {acc} for {species_name}")
        return acc.split('.')[0]  #NCBI FTP path format
    
    except Exception as e:
        print(f"[ERROR] Failed to fetch accession for {species_name}: {e}")
        return None

#download FASTA and GFF from NCBI's FTP:

FTP_HOST = "ftp.ncbi.nlm.nih.gov"
BASE_PATH = "/genomes/all"
output_folder = "user-data"
os.makedirs(output_folder, exist_ok = True) #creates the directory if it doesn't exist already

"""
    This function will fetch FASTA files from NCBI.
    Each FASTA will be saved locally (default: current user director)

    #fetches fasta each time; change to make that optional to user
    
    parameters:
        assemblies (string): NCBI genome assembly accession IDs
        output_folder (string): directory where FASTA files will be stored
"""

def fetch_fasta(species_name, accession, force_download=False):
    accession = accession.split('.')[0]  # remove version (e.g. .1)
    
    # extract only the number part from accession (e.g. '036942975')
    numeric_part = accession.split('_')[1]
    
    # build FTP subdirectory path like 036/942/975
    acc_path = "/".join([numeric_part[i:i+3] for i in range(0, len(numeric_part), 3)])
    
    # final FTP folder path
    full_path = f"{BASE_PATH}/{acc_path}/{accession}"

    print(f"[INFO] Fetching files for {species_name} from {full_path}...")

    try:
        with FTP(FTP_HOST) as ftp:
            ftp.login()
            ftp.cwd(full_path) #will use 'full_path' var from above
            files = ftp.nlst()

            fna_file = next((f for f in files if f.endswith("_genomic.fna.gz")), None) #producing FASTA formatted file
            gff_file = next((f for f in files if f.endswith("_genomic.gff.gz")), None) #producing GFF formatted file 

            for file in [fna_file, gff_file]:
                if file:
                    local_name = f"{species_name.replace(' ', '_')}_{file}"
                    local_path = os.path.join(output_folder, local_name)

                    if os.path.exists(local_path) and not force_download:
                        print(f"[SKIPPED] {file} already exists. Use --force to re-download.")
                        continue

                    with open(local_path, "wb") as f:
                        print(f"[INFO] Downloading {file}...")
                        ftp.retrbinary(f"RETR {file}", f.write)
                    print(f"[SUCCESS] Downloaded {file} to {local_path}")
                else:
                    print(f"[WARNING] File not found for {species_name}")
    except Exception as e:
        print(f"[ERROR] FTP download failed for {species_name}: {e}")

#use argparse for --force option 
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true", help="Force re-download even if files exist")
    return parser.parse_args()

# ==================
# MAIN PROGRAM LOGIC 
# ==================

if __name__ == "__main__":
    args = parse_args()
    species_dict = get_species()

    if not species_dict:
        print("[INFO] No species provided. Exiting.")
    else:
        for species in species_dict.keys():
            accession = get_assembly(species)
            if accession:
                fetch_fasta(species, accession, force_download=args.force)

        print("[COMPLETE] All downloads finished.")






