# pdistmat

A Python script for calculating pairwise nucleotide distances from a multiple sequence alignment in FASTA format.

## Usage

```bash
python pdistmat.py --input alignment.fasta --output distances.tsv
```

The script requires Python 3 and pandas.

## Input

A nucleotide multiple sequence alignment in FASTA format. All sequences must have the same length.

## Output

A tab-separated square distance matrix containing the pairwise genetic distances between all sequences in the alignment.

## Distance calculation

Distances are calculated as the proportion of alignment positions containing incompatible nucleotides.

- Gaps are not counted as differences.
- IUPAC ambiguous nucleotides are considered matches when their possible nucleotide identities overlap.
- All alignment positions, including gaps and ambiguous positions, are retained in the denominator.

## Dependencies

- Python 3
- pandas
