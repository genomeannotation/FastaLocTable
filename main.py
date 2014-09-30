#!/usr/bin/env python
import sys
from src.fasta_reader import read_fasta
from src.translator import translate

def main():
    with open(sys.argv[1], "r") as locations:
        with open(sys.argv[2], "r") as fasta:
            with open(sys.argv[3], "r") as genes_file:
                genes = {}
                repeats = {}
                for line in genes_file:
                    cols = line.strip().split()
                    if len(cols) < 3:
                        continue
                    if cols[2] in genes:
                        if cols[2] in repeats:
                            repeats[cols[2]] += 1
                            cols[2] += "_r"+str(repeats[cols[2]])
                        else:
                            repeats[cols[2]] = 0
                            cols[2] += "_r1"
                    genes[cols[2]] = (int(cols[0]), int(cols[1]), cols[3])
            if len(sys.argv) >= 5:
                seq_extent = int(sys.argv[4])
            else:
                seq_extent = 50
            locations = [int(line) for line in locations]
            seqs = read_fasta(fasta)
            sys.stdout.write("Location\t"+"in_gene\t"+"\t".join([seq.header for seq in seqs])+"\tAminoAcids\tsequence\ttranslation\n")
            sys.stderr.write("Location\t"+"in_gene\t"+"\t".join([seq.header for seq in seqs])+"\tAminoAcids\tsequence\ttranslation\n")
            for location in locations:
                outfile = sys.stdout
                bases = [seq.bases[location-1] for seq in seqs]
                # Condition 1: variation in tail sequences
                if len(set(bases[1:])) < 2:
                    outfile = sys.stderr
                # Condition 2: no more than two different bases
                if len(set(bases)) > 2:
                    outfile = sys.stderr
                # Condition 3: no variation, gaps, or Ns within range of location
                locations_in_range = [x for x in locations if x != location and abs(location-x) <= seq_extent]
                before_seq = seqs[0].bases[location-seq_extent-1:location-1]
                after_seq = seqs[0].bases[location:location+seq_extent]
                if locations_in_range or 'N' in before_seq or 'N' in after_seq or 'n' in before_seq or 'n' in after_seq or '-' in before_seq or '-' in after_seq:
                    outfile = sys.stderr
                unique_bases = list(set(bases))
                # Is SNP in gene?
                in_gene = "."
                strand_map = {">":"+", "<":"-"}
                translation = "."
                amino_acids = "."
                for gene, coords in genes.items():
                    if location >= coords[0] and location <= coords[1]:
                        in_gene = gene
                        snp_gene_rel = location-coords[0] # Coordinate of SNP relative to gene
                        amino_acid_coord = snp_gene_rel // 3 # Coordinate of SNP product amino acid in translated sequence
                        amino_acids = []
                        for seq in seqs:
                            gene_bases = seq.bases[coords[0]-1:coords[1]]
                            t = translate(gene_bases, strand_map[coords[2]])
                            amino_acids.append(t[amino_acid_coord])
                        gene_bases = seqs[0].bases[coords[0]-1:coords[1]]
                        translation = translate(gene_bases, strand_map[coords[2]])
                        translation = translation[:amino_acid_coord]+"["+"/".join(set(amino_acids))+"]"+translation[amino_acid_coord+1:]
                        amino_acids = "/".join(set(amino_acids))
                        break
                outfile.write(str(location)+"\t"+in_gene+"\t"+"\t".join(bases)+"\t"+amino_acids+"\t"+before_seq+"["+unique_bases[0]+"/"+unique_bases[1]+"]"+after_seq+"\t"+translation+"\n")

##################################3

if __name__ == "__main__":
    main()
