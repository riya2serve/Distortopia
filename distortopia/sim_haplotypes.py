#!/usr/bin/env python

"""Simulate haploid genomes variants from a haploid reference genome.

Specify a heterozygosity rate and mutations are randomly uniformly
applied to positions of a haploid genome to mutate to alternate alleles
and return as fasta.

Input: reference.fa
Output: fasta written to stdout
Options: heterozygosity

Usage: $ python sim_haplotypes.py REF HET SEED > hap.fa
"""

import sys
import gzip
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
import numpy as np


def sample_positions(reference, rate, rng):
    """Randomly sample positions between 0 and N at the specified rate
    """
    # get chrom lengths and sumlen    
    chrom_lens = [len(i) for i in reference]
    total_len = sum(chrom_lens)

    # randomly sample positions between 0-N at specified rate
    mask = rng.random(total_len) < rate
    sampled_positions = np.flatnonzero(mask)
    return sampled_positions


def mutate(base, rng):
    """Return a different base than the current one"""
    bases = list(set("ACTG") - set((base)))
    return rng.choice(bases)


def add_variants(reference, heterozygosity, hap_id, seed):
    """Return haplotype with variants mutated relative to reference."""
    rng = np.random.default_rng(seed)

    # assign alt variant to haplotypes
    for chrom in reference:
        varsites = sample_positions(reference[chrom], heterozygosity, rng)
        for site in varsites:
            reference[chrom][site] = mutate(reference[chrom][site], rng)

    # stream fasta to stdout
    records = [SeqRecord(reference[i], id=i, description=f"{hap_id}") for i in sorted(reference)]
    SeqIO.write(records, sys.stdout, "fasta")            


def get_fasta_dict(reference_path):
    """Return a dict mapping seq.id:seq for each chrom in reference
    except for organelles.
    """
    # parse the reference
    if reference_path.name.endswith(".gz"):
        xopen = gzip.open
    else:
        xopen = open
    with xopen(reference_path, "rt") as handle:
        reference = SeqIO.parse(handle, "fasta")
        # convert to a dict[str, MutableSeq]
        reference = {
            i.id: MutableSeq(i.seq.upper()) for i in reference
            if i.id not in ["mitochondria", "chloroplast"]
        }
    return reference



if __name__ == "__main__":

    REF = Path(sys.argv[1])
    HET = float(sys.argv[2])
    SEED = int(sys.argv[3])
    HAP_ID = str(sys.argv[4])

    # parse the reference
    reference = get_fasta_dict(REF)

    # write new haplotype to stdout
    add_variants(reference, HET, HAP_ID, SEED)
