#!/usr/bin/env python
__author__ = 'mahajrod'

import argparse
from RouToolPa.Routines import SequenceRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input",
                    help="Fasta file with sequences")
parser.add_argument("-p", "--prefix", action="store", dest="prefix",
                    help="Prefix of output files")
parser.add_argument("-o", "--output_directory", action="store", dest="output_dir",
                    help="Directory to write output files")
parser.add_argument("-n", "--num_of_records_per_file", action="store", dest="num_of_records_per_file",
                    type=int,
                    help="Number of sequences per output file")
parser.add_argument("-f", "--num_of_out_files", action="store", dest="num_of_out_files",
                    type=int,
                    help="Number of output files")

args = parser.parse_args()

if args.num_of_records_per_file and args.num_of_out_files:
    raise ValueError("Options -n and -f can't be set simultaneously")

SequenceRoutines.split_fasta(args.input, args.output_dir, num_of_recs_per_file=args.num_of_records_per_file,
                             num_of_files=args.num_of_out_files, output_prefix=args.prefix)
