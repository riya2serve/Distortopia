
# Distortopia

Simulating F1 Hybrids and Detecting Segregation Distortion from Genomic Data

## Project Overview
Distortopia is a Python toolkit for simulating haploid F1 hybrid genomes, from diploid parents, and analyzing segregation distortion. Built for evolutionary biologists and genomicists, it automates hybrid genome generation, alignment, long-read simulation, and visual analysis using real SNP positions. 

---

### Current Capabilities
1. Downloads parental ("reference") genomes via NCBI CLI tools
2. Compares parental genomes to detect SNPS and indels
3. Simulates recombinat F1 hybrid genomes with or without crossover events
4. Simulates long reads from F1 genomes and aligns them to parents using `minimap2`
5. Visualizes outputs `.paf`, `.sam`, `.tsv`, and styled `.html` reports


### Installation 
To install the development version locally and contribute to the code:

1. Clone the repository
```bash
git clone https://github.com/riya2serve/Distortopia.git
cd Distortopia
```
You may also install (or clone and build) `minimap2` separately if needed:
```bash
brew install minimap2 #using homebrew (ideal for Mac users)

#OR manually:
git clone https://github.com/lh3/minimap2.git #cloning
cd minimap2
make 
```
2. Install dependencies

With **conda**:
```bash
conda install python=3.10 biopython pandas matplotlib -c conda-forge
conda install -c bioconda minimap2 ncbi-datasets-cli
```
Or with **pip**:
```bash
pip install biopython pandas matplotlib
pip install ncbi-datasets-pylib
```
3. Install the packages in editable mode
```bash
pip install -e .
```

### Usage
Scripts must be run individually from the command line. Example workflows are provided below.

#### 1. Fetch Genomes
Use `choose_fasta()` in `fetch_genomes.py` to download and rehydrate genomes from NCBI

```bash
python populatenome.py --extract
```
- You will be prompted to enter the species names
- Use `--force` to overwrite existing downloads
- Use `--extract` to extract `.fna.gz` and `.gff.gz` files

==== EXAMPLE INTERACTIVE INPUT
```bash
#Species name (or 'done'): 
#Species name (or 'done'):
#...
```
#### 2. Simulate F1 Hybrid
Use `simf1poly.py` to simulate a haploid F1 hybrid genome based on SNP alignment data

```bash
python simf1poly.py \
--ref-dir path/to/parent1.fna \
--query-dir path/to/parent2.fna \
--snp path/to/snp_summary.tsv \
--out F1_hybrid.fna
```
==== EXAMPLE INTERACTIVE INPUT
```bash
#python simf1poly.py \
  #--ref-dir user-data/Arabidopsis_thaliana/genome.fna \
  #--query-dir user-data/Arabidopsis_lyrata/genome.fna \
  #--snp snp_positions.tsv \
  #--out F1_hybrid.fna
```
#### 3. Align Hybrid Reference
Use `aligning_genomes.py` to align simulated hybrid(s) to reference genome(s) using minimap2

```bash
python simf1poly.py \
--ref-dir path/to/parent1.fna \
--query-dir path/to/f1_hybrid.fna \
--out hybrid_vs_parent1.paf \
--mode paf
```
- Default file output is `.paf` (for dotplot visualization)
- You can optionally choose to output `.sam` (for genome browser visualizations; e.g. IGV)

==== EXAMPLE INTERACTIVE INPUT
```bash
#python aligning_genomes.py \
  #--ref-dir user-data/Arabidopsis_thaliana/ncbi_dataset/data/GCA_000001735.2/GCA_000001735.2_TAIR10.1_genomic.fna \
  #--query-dir f1_hybrid.fna \
  #--out hybrid_vs_ref.paf
```

### Future Plans 
- Package Distortopia as a true CLI (python -m distortopia) with `argparse`
- Extend support to simulate multiple haploid F1 hybrids
- Automate SNP calling from minimap2 alignments
- Generate summary plots (recombination patterns and/or rates)
- Improve error handling




