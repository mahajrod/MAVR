#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import SequenceRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-a", "--input_file_A", action="store", dest="input_file_a", required=True,
                    help="Input file A with sequences")
parser.add_argument("-b", "--input_file_B", action="store", dest="input_file_b", required=True,
                    help="Input file A with sequences")
parser.add_argument("-o", "--output_file_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of input and output files. Allowed formats genbank, fasta(default)")

args = parser.parse_args()

SequenceRoutines.compare_sequences_from_files(args.input_file_a, args.input_file_b, args.output_prefix,
                                              format=args.format, verbose=True)



