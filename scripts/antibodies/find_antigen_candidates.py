#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from os import path

from Bio import SeqIO
from Tools.BLAST import BLASTp
from Routines.Sequence import get_kmer_dict_as_seq_records, record_by_expression_generator


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input fasta file with protein sequence")
parser.add_argument("-r", "--protein_id", action="store", dest="protein_id", required=True,
                    help="Id of proteins to ignore in species blast search")
parser.add_argument("-s", "--start", action="store", dest="start", default=1, type=int,
                    help="Start of region of interest in protein(1-based).")
parser.add_argument("-e", "--end", action="store", dest="end", type=int,
                    help="End of region of interest in protein(1-based).")
parser.add_argument("-d", "--db", action="store", dest="species_db", required=True,
                    help="Protein blast db of species")
parser.add_argument("-v", "--species_evalue", action="store", dest="species_evalue", type=float,
                    help="E-value cutoff for species blast search")
parser.add_argument("-u", "--immune_db", action="store", dest="immune_db", required=True,
                    help="Protein blast db of species to be immunized")
parser.add_argument("-m", "--immune_evalue", action="store", dest="immune_evalue", type=float,
                    help="E-value cutoff for immune blast search")
parser.add_argument("-l", "--length", action="store", dest="length", default=15, type=int,
                    help="Length of antigen")
parser.add_argument("-o", "--output_prefix", action="store", dest="out_prefix", default="output",
                    help="Prefix of output files")
parser.add_argument("-t", "--threads", action="store", dest="threads", default=1, type=int,
                    help="Number of threads")
args = parser.parse_args()

kmer_file = "%s.kmer.fasta" % args.out_prefix
species_blast_hits = "%s.species.blast.hits" % args.out_prefix
immune_blast_hits = "%s.immune.blast.hits" % args.out_prefix

BLASTp.threads = args.threads

sequence = list(SeqIO.parse(args.input, format="fasta"))[0]

print("Constructing kmer list...")
#print len(sequence.seq)
kmer_dict = get_kmer_dict_as_seq_records(sequence.seq, args.length, args.start, args.end)
kmer_ids = list(kmer_dict.keys())

SeqIO.write(record_by_expression_generator(kmer_dict), kmer_file, format="fasta")

print("Blast of kmers vs species peptides")
BLASTp.search(kmer_file, args.species_db, outfile=species_blast_hits,
              blast_options=None, evalue=args.species_evalue, output_format=6)
print("Blast of kmers vs immunogenetic species peptides")
BLASTp.search(kmer_file, args.immune_db, outfile=immune_blast_hits,
              blast_options=None, evalue=args.immune_evalue, output_format=6)
