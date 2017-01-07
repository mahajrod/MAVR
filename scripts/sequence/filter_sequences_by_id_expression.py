#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

from Routines import SequenceRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_fasta", action="store", dest="input", required=True,
                    help="Comma separated list of files/directories with sequences")
parser.add_argument("-a", "--filtered_file", action="store", dest="filtered_file", required=True,
                    help="File to write sequences with ids with hits from regular expression")
parser.add_argument("-b", "--filtered_out_file", action="store", dest="filtered_out_file", required=True,
                    help="File to write sequences with ids without hits from regular expression")
parser.add_argument("-e", "--regular_expression", action="store", dest="regular_expression", required=True,
                    help="Regular expression to seek in sequence ids")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of input and output file. Allowed formats genbank, fasta(default)")
parser.add_argument("-p", "--parsing_mode", action="store", dest="parsing_mode", default="index_db",
                    help="Parsing mode for input sequence file. "
                         "Possible variants: 'index_db'(default), 'index', 'parse'")

args = parser.parse_args()

SequenceRoutines.filter_seq_by_reg_expression_from_file(args.input_file, args.regular_expression,
                                                        args.filtered_file, args.filtered_out_file,
                                                        parsing_mode="index_db", format=args.format,
                                                        index_file="tmp.idx", retain_index=False, reg_exp_flags=None)
