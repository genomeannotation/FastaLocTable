#!/usr/bin/env python
import sys
from src.fasta_reader import read_fasta

def main():
    with open(sys.argv[1], "r") as locations, open(sys.argv[2], "r") as fasta:
        locations = [int(line) for line in locations]
        seqs = read_fasta(fasta)
        print("Location\t"+"\t".join([seq.header for seq in seqs]))
        for location in locations:
            print(str(location)+"\t"+"\t".join([seq.bases[location-1] for seq in seqs]))

##################################3

if __name__ == "__main__":
    main()
