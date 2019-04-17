#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse
from Bio import SeqIO
from RouToolPa.Routines.File import make_list_of_path_to_files
from RouToolPa.Routines.Sequence import record_by_expression_generator



parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file_list", action="store", dest="input", required=True,
                    type=lambda s: make_list_of_path_to_files(s.split(",")),
                    help="Comma-separated list of input files/directories with sequences")
parser.add_argument("-o", "--output_file", action="store", dest="output", required=True,
                    help="Output file with sequences")
parser.add_argument("-f", "--input_format", action="store", dest="input_format", required=True,
                    help="Format of input file. Allowed formats genbank, fasta, etc.")
parser.add_argument("-u", "--output_format", action="store", dest="output_format", required=True,
                    help="Format of input file. Allowed formats genbank, fasta, etc.")

args = parser.parse_args()

tmp_index_file = "temp.idx"

print("Parsing %s..." % (args.input if isinstance(args.input, str) else ",".join(args.input)))

sequence_dict = SeqIO.index_db(tmp_index_file, args.input, format=args.input_format)
SeqIO.write(record_by_expression_generator(sequence_dict), args.output, format=args.output_format)
os.remove(tmp_index_file)
