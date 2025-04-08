import os 
from ftplib import FTP
import gzip #for FASTA and GFF file decompression 
from Bio import SeqIO #BioPython for NCBI access and FASTA parsing
import pandas as pd #for dataframe handling (not necessary)

# ================
# CONFIGURATION
# ================

FTP_HOST = "ftp.ncbi.nlm.nih.gov"
BASE_PATH = "/genomes/all"

#Step 1. list genome assembly accessions
assemblies = {
    "A_retroflexus": "JBEFNP000000000", #(data made available by Raieymo et al. 2025)
    "A_hybridus": "JBEL0C000000000" #(data made available by Raieymo et al. 2025)
}

#Step 2. define output directory to store the downloaded FASTA files
output_folder = "user-data"
#creates the directory if it doesn't exist already
os.makedirs(output_folder, exist_ok = True)

# ================
# MY FUNCTIONS
# ================

def fetch_fasta(assembly_name, accession):
    """
    This function will fetch FASTA files from NCBI.
    Each FASTA will be saved locally (default: current user director)
    
    parameters:
        assemblies (string): NCBI genome assembly accession IDs
        output_folder (string): directory where FASTA files will be stored
    """

    acc_path = "/".join([accession[i:i+3] for i in range(0, len(accession), 3)]) 
    #breaking accession into chunks of characters bc that is how NCBI organizes the FTP directory
    full_path = f"{BASE_PATH}/{acc_path}/{accession}"
    #exact FTP folder when NCBI is storing the genome assembly 

    print(f"[INFO] Fetching {assemblies}. . .")

    try:
        with FTP(FTP_HOST) as ftp:
            ftp.login()
            ftp.cwd(full_path) #will use 'full_path' var from above
            files = ftp.nlst()

            fna_file = next((f for f in files if f.endswith("_genomic.fna.gz")), None) #producing FASTA formatted file
            gff_file = next((f for f in files if f.endswith("_genomic.gff.gz")), None) #producing GFF formatted file 

            for file in [fna_file, gff_file]:
                if file:
                    local_path = os.path.join(output_folder, f"{assembly_name}_{file}")
                    with open(local_path, "wb") as f:
                        print(f"[INFO] Downloading {file}...")
                        ftp.retrbinary(f"RETR {file}", f.write)
                    print(f"[SUCCESS] Downloaded {file} to {local_path}")
                else:
                    print(f"[WARNING] Could not find file for {assembly_name}")
    except Exception as e:
        print(f"[ERROR] Failed to download {assembly_name}: {e}")

# ==============
# FETCH CHROMS
# ==============
if __name__ == "__main__":
    for acc_id in accession_ids:
        fetch_fasta(acc_id, output_folder)

    print("[COMPLETE] Finished downloading all chromosomes.")

"""
==================================
EXPECTED OUTPUT (in terminal)
==================================

[INFO] Fetching NC_079487.1 from NCBI...
[SUCCESS] Saved NC_079487.1 to spinach_genome/NC_079487.1.fasta
"""










