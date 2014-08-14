#!/usr/bin/env python
import sys
from src.fasta_reader import read_fasta

def main():
    with open(sys.argv[1], "r") as locations:
        with open(sys.argv[2], "r") as fasta:
            if len(sys.argv) >= 4:
                seq_extent = int(sys.argv[3])
            else:
                seq_extent = 50
            locations = [int(line) for line in locations]
            seqs = read_fasta(fasta)
            print("Location\t"+"\t".join([seq.header for seq in seqs]))
            for location in locations:
                bases = [seq.bases[location-1] for seq in seqs]
                # Condition 1: variation in tail sequences
                if len(set(bases[1:])) < 2:
                    continue
                # Condition 2: no more than two different bases
                if len(set(bases)) > 2:
                    continue
                # Condition 3: no variation, gaps, or Ns within range of location
                locations_in_range = [x for x in locations if x != location and abs(location-x) <= seq_extent]
                before_seq = seqs[0].bases[location-seq_extent-1:location-1]
                after_seq = seqs[0].bases[location:location+seq_extent]
                if locations_in_range or 'N' in before_seq or 'N' in after_seq or 'n' in before_seq or 'n' in after_seq or '-' in before_seq or '-' in after_seq:
                    continue
                unique_bases = list(set(bases))
                print(str(location)+"\t"+"\t".join(bases)+"\t"+before_seq+"["+unique_bases[0]+"/"+unique_bases[1]+"]"+after_seq)

##################################3

if __name__ == "__main__":
    main()
