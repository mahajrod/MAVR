#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
import os

from Bio import SeqIO
from BCBio import GFF

from RouToolPa.Routines import SequenceRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input_file",
                    help="Input file with sequences")
parser.add_argument("-o", "--output_file", action="store", dest="out_file",
                    help="Output gff with gaps")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of input. Allowed formats genbank, fasta(default)")

args = parser.parse_args()

tmp_index_file = "temp.idx"

print("Parsing %s..." % args.input_file)
sequence_dict = SeqIO.index_db(tmp_index_file, args.input_file, format=args.format)
gaps_dict = SequenceRoutines.find_gaps(sequence_dict)

with open(args.out_file, "w") as out_fd:
    GFF.write(SequenceRoutines.record_by_id_generator(gaps_dict, gaps_dict.keys(),
                                                      verbose=True), out_fd)
os.remove(tmp_index_file)
