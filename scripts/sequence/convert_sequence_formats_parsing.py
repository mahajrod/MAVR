#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse
from Bio import SeqIO
from RouToolPa.GeneralRoutines.File import make_list_of_path_to_files
from RouToolPa.Routines.Sequence import record_by_expression_generator



parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input file with sequences")
parser.add_argument("-o", "--output_file", action="store", dest="output", required=True,
                    help="Output file with sequences")
parser.add_argument("-f", "--input_format", action="store", dest="input_format", required=True,
                    help="Format of input file. Allowed formats genbank, fasta, etc.")
parser.add_argument("-u", "--output_format", action="store", dest="output_format", required=True,
                    help="Format of output file. Allowed formats genbank, fasta, etc.")

args = parser.parse_args()

sequence_dict = SeqIO.to_dict(SeqIO.parse(args.input, format=args.input_format))
SeqIO.write(record_by_expression_generator(sequence_dict), args.output, format=args.output_format)
