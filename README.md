
# Distortopia

Distortopia is a command-line interface (CLI) tool to simulate F1 hybrid genomes and detect segregation distortion from diploid parental VCF inputs. It is designed for use with any organism and is useful for geneticists, breeders, and researchers studying hybridization, recombination, or marker bias.

---

### In development

To install the development version locally and contribute to the code:

1. Clone the repository
```bash
git clone https://github.com/YOUR-USERNAME/Distortopia.git
cd Distortopia
```

2. Install dependencies (user chooses one)

With **conda**:
```bash
conda install python=3.10 biopython pandas matplotlib -c conda-forge
```
with **pip**:
```bash
pip install biopython pandas matplotlib
pip install ncbi-datasets-pylib
```

3. Install the package in editable mode
```bash
pip install -e .
```

4. Run the CLI tool:
```bash
python -m distortopia --help
python populatenome.py --extract
unzip user-data/Arabidopsis_lyrata.zip -d user-data/Arabidopsis_lyrata
unzip user-data/Arabidopsis_thaliana.zip -d user-data/Arabidopsis_thaliana
```
```bash
datasets rehydrate --directory user-data/Arabidopsis_lyrata
datasets rehydrate --directory user-data/Arabidopsis_thaliana
```

For detailed goals, usage examples, pseudocode, and file formats, please see: [proposal.md](./proposal.md).




