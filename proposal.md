
# Distortopia: 
Simulating F1 Hybrids and Detecting Segregation Distortion from Genomic Data

## What task/goal does this project accomplish and why is this useful?

Distortopia provides a modular command-line toolkit for simulating F1 hybrid genomes and detecting segregation distortion using real or simulated genomic data. Segregation distortion (where allele frequencies deviate from Mendelian expecatations) can indicate the presence of (a)selection, (b) hybrid incompatibility, or (c) meitoic drive. These possibilities make this toolkit valuable for evolutionary biologists, geneticists, and breeders. 

Unlike existing tools, Distortopia integrates hybrid genome simulation, sequence alignment, and distortion analysis into a single, customizable pipeline. It supports general-purpose use on any organism with available genome assemblies in FASTA format. 

## What type of data/input should users provide to the program?

Users can input:
- Two parent genome assemblies in FASTA format
- SNP position tables (generated or supplied by the user as a `.tsv`)
- *(Optional)* Pre-existing hybrid genome assemblies or simulated VCFs

These files can be downloaded using built-in modules or prepared independently

## Where will the user data come from?

Distortopia includes scripts to download and prepare reference genome assemblies directly from NCBI. Data can also be sourced from the user's existing genomic datasets (e.g. VCF, experimental genomes).

Example file formats:

```**VCF**
##fileformat=VCFv4.2
##source=Distortopia
#CHROM  POS    ID      REF     ALT     QUAL   FILTER  INFO FORMAT parent1
chr1    10523   .       G       A       .      PASS    .    GT    0/1
chr1    20831   .       T       C       .      PASS    .    GT    1/1
chr1    31005   .       A       G       .      PASS    .    GT    0/0
```
or as:
```**FASTA**
>chr1 length=30427671 assembly=GCF_000001405.39 source=Homo_sapiens
AGCTGACCTAGGCTACCTTACGATCGATCGATCGATCGATGCTAGCTAGCTAGCTGATCGATCGATCGATCGA
CGATCGATCGTACGTACGTAGCTAGCTAGCTAGCTAGCATCGATCGATCGATCGATCGTAGCTAGCTAGCTAG

>chr2 length=242193529 assembly=GCF_000001405.39 source=Homo_sapiens
TTGACGATCGATCGTACTGACTGATCGATCGATAGCTAGCTAGCTAGCTAACGTAGCTAGCTAGCTAGCTAGC
GATCGATCGTAGCTAGCTGATCGATCGATGATCGTAGCTAGCTAGCATCGATCGATCGATCGATCGATCGTAG
```
## How will users interact with the program?

This package is CLI-driven. Modules can be executed individually or through a unified command structure via `__main__.py.` 

Example usage includes:

#Generate an F1 Hybrid
```**bash**
python distortopia/simf1poly.py \
  --ref-dir distortopia/user-data/Arabidopsis_thaliana/...fna \
  --query-dir distortopia/user-data/Arabidopsis_lyrata/...fna \
  --snp distortopia/genomes/snp_positions.tsv \
  --out distortopia/genomes/f1_hybrid.fna
```
#Align F1 genome to reference
```**bash**
python distortopia/aligning_genomes.py \
  --ref-dir distortopia/user-data/Arabidopsis_thaliana/...fna \
  --query-dir distortopia/genomes/f1_hybrid.fna \
  --out distortopia/genomes/hybrid_vs_ref.paf
```

## What type of output will the program produce? 

Distortopia generates:

- Simulated hybrid genome FASTA files
- `.paf` alignments between genomes (via minimap2)
- SNP position TSV files for hybrid construction
- Interactive HTML summaries of SNPs (`athal_vs_alyr_summary.html`)
- Diagnostic dotplot-style PNG visualizations of genome alignments (via D-Genies)

These outputs are suitable for validation, visualization, or identifying distorted genomic regions.

## What other tools currently exist to do this task, or something similar?

While a few tools exist for simulating recombination, generating hybrid genomes, or analyzing variant data, none fully integrates simulation and segregation distortion into a lightweight, CLI-friendly package. Below is a comparative overview: 

[recom-sim](https://github.com/salanova-elliott/recom-sim):
A Python-based simulator that generates recombinant genomes using parental genotypes and genetic maps. It excels in modeling crossover interference and recombination rate variation. However, it does not support distortion detection, marker parsing, or hybrid genome output without a genetic map. This limits its generality compared to Distortopia.

[simuPOP](https://github.com/BoPeng/simuPOP):
A powerful forward-time simulator designed for modeling complex demographic and evolutionary dynamics. While it is extremely flexible, its steep learning curve and scripting-heavy interface make it less accessible for quick hybrid genome simulations and/or small-scale empirical testing.

[vcftools](https://github.com/vcftools/vcftools) and [bcftools](https://github.com/samtools/bcftools):
These standard bioinformatics tools are essential for VCF filtering, annotation, and manipulation. However, they do not provide any genome simulation or segregation distortion functionality. They are best used as pre-processing or post-processing utilities in workflows like Distortopia's 

# What make Distortopia unique?

- Simulates haploid F1 hybrid genomes directly from base-level alignments and SNP data
- Detects and exports SNP positions with potential segregation distortion
- Uses minimap2 `.paf` alignments as a foundation for simulation, not variant calls or maps
- Supports both real and simulated data inputs
- Outputs are ready for visualization (HTML summaries, dotplots, alignments) or use in downstream analyses
- Designed for ease of use

Ultimately, Distortopia fills a current gap in population genomic tools: the ability to generate and analyze hybrid genomes **without needing**  phased VCFs, recombination maps, or full pedigree simulations. This makes it well-suited for researchers working with non-model organissm, educators teaching hybridization concepts, or developers testing alignment and distortion pipelines.

### Psuedocode 

Each component of Distortopia is modularized in its own Python script. Each script contains detailed docstrings, inline comments, and function-level documentation that serves as both pseduocode and developer guidance. These scropts form the backbone of the project:

- `populatenome.py`:
	Interactively selects FASTA genomes from two species, and aligns them 
- `simparents.py`:

- `simf1poly.py`:

- `aligning_genomes.py`:
	Aligns a hybrid genome back to a parental reference using `minimap2`, optionally in `.paf` or `.sam` format. Useful for validating F1 simulation and visualizing genome-wide similarity (e.g. via dot plots on [D-genesis](https://dgenies.toulouse.inra.fr).
- `module.py`: (*optional utility*)
	Houses helper functions
- `__main__.py`:
	Intended for future CLI integration. This will act as a wrapper allowing users to invoke the full pipeline with flags like `--simulate-f1` or `--detect-distortion`.


