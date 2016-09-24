#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse

from Routines.File import make_list_of_path_to_files
from Routines import SequenceRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file_list", action="store", dest="input", required=True,
                    type=lambda s: make_list_of_path_to_files(s.split(",")),
                    help="Comma-separated list of input files/directories with sequences")
parser.add_argument("-o", "--output_file", action="store", dest="output", default="stdout",
                    help="Output file with sequences")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of input and output files. Allowed formats genbank, fasta(default)")
parser.add_argument("-d", "--id_file", action="store", dest="id_file",
                    help="File with ids of sequences to extract")
parser.add_argument("-c", "--id_column", action="store", dest="id_column", type=int, default=0,
                    help="Number of column with ids in id file (0-based). Default: 0")

args = parser.parse_args()

SequenceRoutines.extract_sequence_by_ids(args.input, args.id_file, args.output, format=args.format, verbose=True,
                                         id_column_number=args.id_column)
