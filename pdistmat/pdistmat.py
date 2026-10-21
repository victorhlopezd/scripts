#!/usr/bin/env python3

import argparse
import pandas as pd

def fasta_to_dict(file):
    seqs = {}
    with open(file) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:]
                seqs[name] = ""
            else:
                seqs[name] += line.upper()
    return seqs


def distance(seq1, seq2):
    ambigs = {
        "A": "A", "C": "C", "G": "G", "T": "T",
        "R": "AG", "Y": "CT", "S": "GC", "W": "AT",
        "K": "GT", "M": "AC", "B": "CGT", "D": "AGT",
        "H": "ACT", "V": "ACG", "N": "ACGT"
    }

    dist = 0

    for n1, n2 in zip(seq1, seq2):
        if n1 == "-" or n2 == "-":
            continue
        if not set(ambigs[n1]).intersection(ambigs[n2]):
            dist += 1

    return dist / len(seq1)


def distmat(file, output):
    alignment = fasta_to_dict(file)

    matrix = pd.DataFrame(index=alignment.keys(), columns=alignment.keys())

    for seq1 in alignment:
        for seq2 in alignment:
            matrix.loc[seq1, seq2] = distance(alignment[seq1], alignment[seq2])

    matrix.to_csv(output, sep="\t")


parser = argparse.ArgumentParser(
    description="Calculate pairwise nucleotide distances from an aligned FASTA file. "
                "Gaps are ignored and ambiguous nucleotides are considered matches "
                "when their possible nucleotide identities overlap."
)

parser.add_argument(
    "--input", "-i",
    required=True,
    help="Input multiple sequence alignment in FASTA format"
)

parser.add_argument(
    "--output", "-o",
    required=True,
    help="Output TSV distance matrix"
)

args = parser.parse_args()

distmat(args.input, args.output)