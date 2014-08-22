#!/usr/bin/env python

import sys
from src.fasta_reader import read_fasta

def to_gapped_index(gaps, index):
    for gap in gaps:
        if gap <= index:
            index += 1
        else:
            return index
    return index

def main():
    with open(sys.argv[1], "r") as coords:
        with open(sys.argv[2], "r") as gapped_fasta:
            gapped_seqs = read_fasta(gapped_fasta)
            # Gaps = [index]
            gaps = []
            bases = gapped_seqs[3].bases
            for index, base in enumerate(bases):
                if base == "-":
                    gaps.append(index+1)
            # Genes = [(start, stop, gene_name, strand)]
            genes = []
            for line in coords:
                cols = line.strip().split()
                # Skip if it's not a gene line
                if cols[0] != ">" and cols [0] != "<":
                    continue
                start = int(cols[1])
                end = int(cols[2])
                name = cols[3]
                strand = cols[0]
                # Convert to gapped space
                start = to_gapped_index(gaps, start)
                end = to_gapped_index(gaps, end)
                # Add gene
                genes.append((start, end, name, strand))
            with open("genes.tsv", "w") as genes_file:
                for gene in genes:
                    genes_file.write("\t".join([str(g) for g in gene])+"\n")

###############################

if __name__ == "__main__":
    main()
