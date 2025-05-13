
# Distortopia

Simulating F1 Hybrids and Detecting Segregation Distortion from Genomic Data

## Project Overview
Distortopia is a Python toolkit for simulating haploid F1 hybrid genomes, from diploid parents, and analyzing segregation distortion. Built for evolutionary biologists and genomicists, it automates hybrid genome generation, alignment, long-read simulation, and visual analysis using real SNP positions. 

---

### Current Capabilities
1. Downloads parental ("reference") genomes via NCBI CLI tools
2. Compares parental genomes to detect SNPS and indels
3. Simulates recombinat F1 hybrid genomes with or without crossover events
4. Simulates long reads from F1 genomes and aligns them to both parents using `minimap2`
5. Parses alignments into VCF-style TSVs and HTML summaries
6. Visualizes output reports

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

#### 1. Download Parental Genomes
Use `scripts/fetch_genomes.py`

```bash
python scripts/fetch_genomes.py --extract
```
- You will be prompted to enter the species names
- Use `--force` to overwrite existing downloads
- Use `--extract` to extract `.fna.gz` and `.gff.gz` files

==== EXAMPLE INTERACTIVE INPUT
```bash
#Species name (or 'done'): Arabidopsis thaliana
#Species name (or 'done'): Arabidopsis lyrata
#Species name (or 'done'): done
#...
```
#### 2. Compare Parent Genomes and Generate SNP summary
Use `scripts/align_parents.py `

```bash
python simf1poly.py \
--ref-dir path/to/parent1.fna \
--query-dir path/to/parent2.fna \
--snp path/to/snp_summary.tsv \
--out F1_hybrid.fna
```
==== EXAMPLE INTERACTIVE INPUT
```bash
python scripts/align_parents.py \
  --ref-dir input-data/Arabidopsis_thaliana \
  --alt-dir input-data/Arabidopsis_lyrata \
  --outdir genomes/par_alignments \
  --threads 8
```
This outputs:
- `.paf`: Minimap2 alignment summary
- `.html`: Visual summary
- `.snp_positions.tsv`: SNPs coordinates (necessary for F1 simulation)
```bash
open genomes/athal_vs_alyr_summary.html
```
#### 3. Simulate Recombinant F1 hybrids 
Use `scripts/create_f1_hybrid.py`

```bash
python scripts/create_f1_hybrid.py \
  --ref-fasta path/to/ref.fna \
  --alt-fasta path/to/alt.fna \
  --snp-tsv genomes/par_alignments/snp_positions.tsv \
  --outdir genomes/f1_simulations \
  --reps 3
```
This outputs: 
- `f1_genome_repx.fna`, `f1_table.csv`
To reduce scope for debugging and avoiding error:
```bash
head -n 10 genomes/par_alignments/snp_positions.tsv > genomes/par_alignments/snp_subset.tsv
```
#### 4. Simulate Long Reads 
Use `scripts/long_reads.py`

```bash
python scripts/long_reads.py \
  --ref-fasta input-data/Arabidopsis_thaliana/ncbi_dataset/data/GCA_020911765.2/GCA_020911765.2_ASM2091176v2_genomic.fna \
  --alt-fasta input-data/Arabidopsis_lyrata/ncbi_dataset/data/GCA_000524985.1/GCA_000524985.1_Alyr_1.0_genomic.fna \
  --snp-tsv genomes/par_alignments/snp_positions.tsv \
  --outdir genomes/f1_simulations \
  --reps 1
```
This outputs:
- `f1_reads_repx.fasta`
Want to know the total number of long reads generated? Try this:
```bash
grep ">" genomes/f1_simulations/f1_reads_rep1.fasta | wc -l
```
==== EXAMPLE INTERACTIVE INPUT
```bash
#Field            What it tells you                                   Source
##POS             Read alignment start on reference genome            SAM field 4
##crossover_pos   Breakpoint between parental segments in F1 genome   f1_table.csv
```
#### 5. Map Long Reads to Reference Genomes
Use `scripts/map_long_reads.py`

```bash
python scripts/map_long_reads.py \
  --ref-fasta input-data/Arabidopsis_thaliana/ncbi_dataset/data/GCA_020911765.2/GCA_020911765.2_ASM2091176v2_genomic.fna \
  --alt-fasta input-data/Arabidopsis_lyrata/ncbi_dataset/data/GCA_000524985.1/GCA_000524985.1_Alyr_1.0_genomic.fna \
  --snp-tsv genomes/par_alignments/snp_positions.tsv \
  --outdir genomes/f1_simulations \
  --reps 1 \
  --threads 8
```
This outputs:
- `repx_vs_alt.sam`, `repx_vs_ref.sam`

#### 6. Parse `.sam` outputs to a VCF-like Output Summary File
Use `scripts/sam_to_vcf.py`

```bash
python scripts/sam_to_vcf.py \
  --sam-ref distortopia/genomes/f1_simulations/rep1_vs_ref.sam \
  --sam-alt distortopia/genomes/f1_simulations/rep1_vs_alt.sam \
  --output distortopia/genomes/f1_simulations/rep1_long_read_alignment.vcf.tsv \
  --rep 1
```
This produces a TSV output that summarizes real SNP positions, REF/ALT alleles, and alignment metadata

#### 7. Visualize F1 Mapping Output
Use `scripts/visualize_f1_csv.py`

```bash
python scripts/visualize_f1_csv.py \
  --input genomes/f1_simulations/rep1_ref_variants.tsv \
  --f1-table genomes/f1_simulations/f1_table.csv
```

### Future Plans 
- Package Distortopia as a true CLI (python -m distortopia) with `argparse`
- Extend support to simulate multiple haploid F1 hybrids
- Automate SNP calling from minimap2 alignments
- Generate summary plots (recombination patterns and/or rates)
- Improve error handling




