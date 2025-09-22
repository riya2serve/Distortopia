#!/usr/bin/env python

"""Simulate long reads from a pool of synthetic gametes.

Given two haploid genomes recombinant gametes are generated from which
long reads are simulated. To accommodate divergent genomes, we use
minimap2 to find homologous regions between genomes where recombination
can occur under the assumption that recombination does not occur at
unpaired regions that have poor homology.

TODO
----
- The trick here seems to be to create a set of pairing rules that should
be based on biology about how we think chromosomes pair up. In other
words, homology at repetitive regions may be just as important or more
than homology at more unique regions. 
- Maybe we should first reduce the pairwise alignment to blocks of 
realistic pairs.
- change mode: sample a position on reference, then find its best match


Input: REF1.fa REF2.fa
Output: fastq to stdout

Usage: $ python sim_gamete_pool_longreads.py REF1 REF2 > gametes.fq
"""

from typing import List, Tuple
import sys
import gzip
import numpy as np
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
import numpy as np


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
        reference = {
            i.id: i.seq.upper() for i in reference
            if i.id not in ["mitochondria", "chloroplast"]
        }
    return reference


def sample_homologous_positions(
    paf_path,
    min_identity=0.95,
    min_aln_len=1000,
    min_s1_s2_ratio=1.5,
    nreps=1,
    rng=None
):
    """Returns a list of tuples of homologous positions
    Args:
        min_s1_s2_ratio: float, minimum ratio s1 / s2 to retain a record
    """
    if rng is None:
        rng = np.random.default_rng()

    homologs_by_ref = defaultdict(list)

    with open(paf_path) as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue

            # Core fields
            ref_name = fields[0]
            ref_len = int(fields[1])
            ref_start = int(fields[2])
            ref_end = int(fields[3])
            strand = fields[4]
            query_name = fields[5]
            query_len = int(fields[6])
            query_start = int(fields[7])
            query_end = int(fields[8])
            match_len = int(fields[9])
            aln_len = int(fields[10])

            # Optional fields (PAF tags)
            tag_dict = {}
            for tag in fields[12:]:
                try:
                    key, typecode, val = tag.split(":", 2)
                    tag_dict[key] = val
                except ValueError:
                    continue

            s1 = int(tag_dict.get("s1", match_len))
            s2 = int(tag_dict.get("s2", "0"))

            # Check identity and s1/s2 ratio
            identity = s1 / aln_len
            if (
                identity >= min_identity and
                aln_len >= min_aln_len and
                ref_end > ref_start and
                query_end > query_start and
                (s2 == 0 or s1 >= min_s1_s2_ratio * s2)
            ):
                homologs_by_ref[ref_name].append({
                    "ref_start": ref_start,
                    "ref_end": ref_end,
                    "query_start": query_start,
                    "query_end": query_end,
                    "query_name": query_name,
                    "strand": strand
                })

    # TODO: return a generator
    results = []
    for _ in range(nreps):
        rep_result = []
        for ref_name, regions in homologs_by_ref.items():
            if not regions:
                continue
            region = rng.choice(regions)
            ref_pos = rng.integers(region["ref_start"], region["ref_end"])
            ref_offset = ref_pos - region["ref_start"]

            if region["strand"] == "+":
                query_pos = region["query_start"] + ref_offset
            else:
                query_pos = region["query_end"] - ref_offset - 1

            rep_result.append((ref_name, ref_pos, region["query_name"], query_pos, region["strand"]))
        results.append(rep_result)

    return results


def iter_simulate_recombinant_genome(ref_seqs, query_seqs, crossover_points, read_len, rng):
    """
    """
    recombinant_records = []

    # iterate over random sampled homologous positions between genomes
    uid = 0
    for i, homologs in enumerate(crossover_points):
        for chrom_match in homologs:
            ref_chr, ref_pos, query_chr, query_pos, strand = chrom_match

            # check that scaffold names are valid
            if ref_chr not in ref_seqs or query_chr not in query_seqs:
                print(f"Warning: {ref_chr} or {query_chr} not found in input FASTAs.")
                continue

            # get the scaffold/chrom from each ref
            ref_seq = ref_seqs[ref_chr]
            query_seq = query_seqs[query_chr]

            # get the sequences from either side of crossover
            if rng.binomial(1, 0.5):
                left = ref_seq[:ref_pos]
                right = query_seq[query_pos:]
            else:    
                left = query_seq[:query_pos]
                right = ref_seq[ref_pos:]

            # reverse complement if strand reverse matched
            if strand == '-':
                right = right.reverse_complement()

            # create the new seq
            new_seq = left + right

            # randomly sample one 'sequenced long read' from this sequence
            start = rng.integers(read_len, len(new_seq) - read_len)
            new_seq = new_seq[start:start + read_len]
            if rng.binomial(1, 0.5):
                new_seq = new_seq.reverse_complement()            

            # is the read recombinant
            recomb = (ref_pos > start) & (ref_pos < start + read_len)

            # label as 'non-recombinant', 'recombinant-unsampled', or 'recombinant-sampled'
            new_id = f"u{uid}_g{i}_recomb{int(recomb)}_ref={ref_chr}:{ref_pos}_q={query_chr}:{query_pos}"
            uid += 1

            # store for writing
            record = SeqRecord(new_seq, id=new_id, description="")
            # record.letter_annotations["phred_quality"] = [40] * len(record)
            recombinant_records.append(record)

        if not i % 1000:
            print(f"simulated {i} genomes")
            yield recombinant_records
            recombinant_records = []
    yield recombinant_records


def simulate_recombinant_genome(ref_seqs, query_seqs, crossover_points, read_len, outfile, rng):
    """Write to fasta in chunks from sim iterator"""
    with open(outfile, 'w') as out:
        for chunk in iter_simulate_recombinant_genome(ref_seqs, query_seqs, crossover_points, read_len, rng):
            SeqIO.write(chunk, out, "fasta")
            out.flush()



def parse_command_line():
    REF1 = Path(sys.argv[1])
    REF2 = Path(sys.argv[2])
    PAF = Path(sys.argv[3])
    SEED = int(sys.argv[4])
    return REF1, REF2, PAF, SEED    


if __name__ == "__main__":

    REF1 = Path("../data/TAIR9_chr_all.fas")
    REF2 = Path("../data/TAIR9_chr_all.fas")    
    # REF2 = Path("../data/Alyrata_384_v1.fa.gz")
    PAF = Path("../data/ref11xref12.paf")
    SEED = 123
    READ_LEN = 20_000
    NGAMETES = int(1e5)

    # REF1, REF2, PAF, SEED = parse_command_line()

    # parse fastas
    ref1 = get_fasta_dict(REF1)
    ref2 = get_fasta_dict(REF2)    

    TOTAL_BASES = (READ_LEN * NGAMETES)
    GENOME_SIZE = sum(len(i) for i in ref1.values())
    GENOME_COV = (TOTAL_BASES / GENOME_SIZE) / 2.
    print(
        f"simulating {NGAMETES:.0f} gametes w/ {READ_LEN/1e4:.0f} Kb "
        f"long-reads for a genome of size {GENOME_SIZE/1e9:.2f}Gb "
        f"to achieve {GENOME_COV:.2f}X coverage")
    # raise SystemExit(0)

    # sample positions
    rng = np.random.default_rng(SEED)
    positions = sample_homologous_positions(
        PAF, 
        min_identity=0.95,
        min_aln_len=20000,
        min_s1_s2_ratio=5,
        nreps=NGAMETES,
        rng=rng,
    )
    # for pos in positions:
    #     print(pos)

    simulate_recombinant_genome(
        ref1, ref2, positions, 
        read_len=READ_LEN, 
        outfile="../data/gametes.fa",
        rng=rng,
    )
    print("done")

    # simulate recomb gametes
    # pos = 46983
    # left1 = ref1['1'][pos-100:pos]
    # left2 = ref2['1'][pos-100:pos]
    # right1 = ref1['1'][pos:pos+100]
    # right2 = ref2['1'][pos:pos+100]

    # print(f"{left1}^{right1}")
    # print(f"{left1}^{right2}")