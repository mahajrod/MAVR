#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse
from Bio import SeqIO
from RouToolPa.Tools.BLAST import BLASTp
from RouToolPa.Routines.Sequence import get_kmer_dict_as_seq_records, record_by_expression_generator




parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input fasta file with protein sequence")
parser.add_argument("-r", "--protein_ids", action="store", dest="protein_ids", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of ids proteins to ignore in species blast search")
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
species_blast_hits_no_self_hits = "%s.species.no_self_hits.blast.hits" % args.out_prefix
species_blast_hits_no_self_hits_ids = "%s.species.no_self_hits.blast.ids" % args.out_prefix

immune_blast_hits = "%s.immune.blast.hits" % args.out_prefix
immune_blast_hits_ids = "%s.immune.blast.ids" % args.out_prefix

nonspecific_fragments_ids = "%s.fragments_with_hits.ids" % args.out_prefix
bad_antigen_candidates_coordinates = "%s.bad_candidates.coordinates" % args.out_prefix
bad_antigen_candidates_coordinates_sorted = "%s.bad_candidates.sorted.coordinates" % args.out_prefix

BLASTp.threads = args.threads

sequence = list(SeqIO.parse(args.input, format="fasta"))[0]

print("Constructing kmer list...\n")
#print len(sequence.seq)
kmer_dict = get_kmer_dict_as_seq_records(sequence.seq, args.length, args.start, args.end)
kmer_ids = list(kmer_dict.keys())

SeqIO.write(record_by_expression_generator(kmer_dict), kmer_file, format="fasta")

print("Blast of kmers vs species peptides\n")
BLASTp.search(kmer_file, args.species_db, outfile=species_blast_hits,
              blast_options=None, evalue=args.species_evalue, output_format=6)

species_grep_string = "grep -v %s %s > %s" % ("|".join(args.protein_ids), species_blast_hits,
                                              species_blast_hits_no_self_hits)
species_awk_string = "awk '{print$1}' %s | uniq > %s" % (species_blast_hits_no_self_hits,
                                                         species_blast_hits_no_self_hits_ids)

os.system(species_grep_string)
os.system(species_awk_string)

print("Blast of kmers vs immunogenetic species peptides\n")
BLASTp.search(kmer_file, args.immune_db, outfile=immune_blast_hits,
              blast_options=None, evalue=args.immune_evalue, output_format=6)
immune_awk_string = "awk '{print$1}' %s | uniq > %s" % (immune_blast_hits,
                                                        immune_blast_hits_ids)
os.system(immune_awk_string)

cat_string = "cat %s %s > %s" % (species_blast_hits_no_self_hits_ids, immune_blast_hits_ids, nonspecific_fragments_ids)
os.system(cat_string)

with open(nonspecific_fragments_ids, "r") as in_fd:
    with open(bad_antigen_candidates_coordinates, "w") as out_fd:
        out_fd.write("#start\tend\n")
        for line in in_fd:
            start, stop = line.strip().split("_")[-1].split("-")
            out_fd.write("%s\t%s\n" % (start, stop))

sort_string = "awk 'NR == 1; NR > 1 {print $0 | \"sort -k1n -k2n\"}'  %s > %s" % (bad_antigen_candidates_coordinates, bad_antigen_candidates_coordinates_sorted)
os.system(sort_string)




