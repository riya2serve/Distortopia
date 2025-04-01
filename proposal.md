
# Distortopia: 
Simulating F1 Hybrids and Detecting Segregation Distortion from Genomic Data

## What task/goal does this project accomplish and why is this useful?

This project provides users with a command-line tool to simulate F1 hybrid genomes and detect genomic regions with segregation distortion. Segregation distortion can indicate the presence of (a) selection, (b) incompatibility, or (c) meiotic drive making this a useful tool for population geneticists, breeders, and anyone studying hybridization.

The project is designed to be general-purpose: any organism with variant data in VCF format and reference sequences in FASTA format can be used. This makes the tool widely applicable for simulation-based studies or exploring real hybrid genotypes.

## What type of data/input should users provide to the program?

The user will provide:
- Two parent genotype files in VCF format (`parent1.vcf`, `parent2.vcf`)
- Reference genomes or chromosome assemblies in FASTA format
- *(Optional)* A real (or simulated) F1 VCF file from which to directly analyze segregation distortion

Command-line arguments will allow users to toggle with the simulation and analysis steps.

## Where will the user data come from?

Data files can be downloaded from public repositories (NCBI, Ensembl) or generated using variant calling pipelines such as GATK, bcftools, or samtools.  

Example file formats:

```**VCF**:
##fileformat=VCFv4.2
##source=Distortopia
#CHROM  POS    ID      REF     ALT     QUAL   FILTER  INFO FORMAT parent1
chr1    10523   .       G       A       .      PASS    .    GT    0/1
chr1    20831   .       T       C       .      PASS    .    GT    1/1
chr1    31005   .       A       G       .      PASS    .    GT    0/0
```
```**FASTA**:
>chr1 length=30427671 assembly=GCF_000001405.39 source=Homo_sapiens
AGCTGACCTAGGCTACCTTACGATCGATCGATCGATCGATGCTAGCTAGCTAGCTGATCGATCGATCGATCGA
CGATCGATCGTACGTACGTAGCTAGCTAGCTAGCTAGCATCGATCGATCGATCGATCGTAGCTAGCTAGCTAG

>chr2 length=242193529 assembly=GCF_000001405.39 source=Homo_sapiens
TTGACGATCGATCGTACTGACTGATCGATCGATAGCTAGCTAGCTAGCTAACGTAGCTAGCTAGCTAGCTAGC
GATCGATCGTAGCTAGCTGATCGATCGATGATCGTAGCTAGCTAGCATCGATCGATCGATCGATCGATCGTAG
```

## How will users interact with the program?

Users will interact with the program via a command-line interface (CLI). To install and run the program locally, they can follow these steps:

1. **Clone the Repository**

Navigate to the current working directory where you intent to clone the project, then run:
```bash
git clone https://github.com/YOUR-USERNAME/Distortopia.git
cd Distortopia
```

2. **Install Dependencies**

If using **conda**:
```bash
conda install python=3.10 biopython -c conda-forge
```
Alternatively, via **pip**:
```bash
pip install biopython
```
Additionaly dependencies (pandas, matplotlib) may be required depending on which modules are used.

3. **Run the program:**

 Execute the program using the -m flag or directly call on specific scripts. 

 For example:
```bash
python -m distortopia --simulate-f1 --parent1 data/parent1.vcf --parent2 data/parent2.vcf --output data/f1_hybrid.vcf

python -m distortopia --detect-distortion --f1 data/f1_hybrid.vcf --output results/segdist_table.csv
```
```bash
optional arguments:
 -h, --help                 show this help message and exit
 --parent1 PARENT1          path to first parental file (FASTA/VCF)
 --parent2 PARENT2          path to second parental file (FASTA/VCF)
 --output OUTPUT            path to the simulated F1 chromosome VCF file
 --snp-count SNP_COUNT  number of SNPs to simulate per chromosome (default: 1000)
 --chrom-length CHROM_LENGTH  chromosome length in base pairs (default: 1Mb)
 --distortion                      apply segregation distortion 
```

## What type of output will the program produce? 

The program will be able to generate:

- Simulated VCF files of F1 hybrids
- CSV tables showing observed vs. expected genotype frequencies
- Diagnostic plots (histograms, barplots) for distorted regions
- *(Optional)* Summary stats exported to JSON or HTML

These outputs are suitable for visualization, downstream analysis, or direct reporting in scientific works.

## What other tools currently exist to do this task, or something similar?

Few tools exist for simulating recombination, hybrid genomes, or exploring genetic variation patterns. Notable examples include:

[recom-sim](https://github.com/salanova-elliott/recom-sim):
This Python-based simulator creates recombinant genomes using parental input genotypes and genetic maps. It is specially designed for modeling recombination under different biological scenarios, including crossover interference. However, the program focuses on recombination mechanics and does not include modules for segregation distortion detection or downstream variant filtering and visualization. Additionally, it assumes the use of a recombination map, whereas Distortopia simulates hybrids and distortion from VCF inputs alone, making it more general-purpose.

[simuPOP](https://github.com/BoPeng/simuPOP):
A forward-time simulator capable of modeling complex demographic and evolutionary scenarios. While the program is extremely flexible, its steep learning curve and scripting-heavy interface make it less accessible for quick hybrid simulations.

[vcftools](https://github.com/vcftools/vcftools) and [bcftools](https://github.com/samtools/bcftools):
These are standard tools for filtering, summarizing, and manipulating VCF files. While they are useful for general variant processing, they do not offer simulation or distortion-specific workflows.

Distortopia is unique in that it can:

- Use real or simulated VCF input from any diploid organism
- Simulate hybrid genotypes
- Detect markers with segregation distortion
- Produce both summary tables and visualizations

The program is intended to be lightweight, customizable, and CLI-driven with minimal user setup. To the author's knowledge, no current tool offers *this* combination of features, making it an extremely useful utility for geneticists, educators, and bioinformatics researchers alike.

### Psuedocode 

Pseudocode has been written directly into the corresponding Python file modules within the distortopia/ directory. Each script includes docstrings and inline comments that outline the purpose of its functions or classes. This includes:

- __main__.py for command-line interface handling

- __simf1poly__.py for F1 hybrid simulation from VCFs

- __segdistorters__.py for detecting segregation distortion

- [...]

These serve as the backbone of the project and will be iteratively converted into fully functional code. This approach will ensure modular development and allow others to immediately begin contributing to and/or testing specific components.

