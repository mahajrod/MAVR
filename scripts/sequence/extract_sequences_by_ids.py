#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
import os

from Bio import SeqIO

from Routines.File import read_ids
from Routines.Sequence import record_by_id_generator


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input_file",
                    help="Input file with sequences")
parser.add_argument("-o", "--output_file", action="store", dest="output_file",
                    help="Output file with sequences")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of input and output files. Allowed formats genbank, fasta(default)")
parser.add_argument("-d", "--id_file", action="store", dest="id_file",
                    help="File with ids of sequences to extract")

args = parser.parse_args()

tmp_index_file = "temp.idx"

id_list = read_ids(args.id_file)

print("Parsing %s..." % args.input_file)
#print("Not found:")
sequence_dict = SeqIO.index_db(tmp_index_file, args.input_file, format=args.format)
SeqIO.write(record_by_id_generator(sequence_dict, id_list), args.output_file, format=args.format)
os.remove(tmp_index_file)



