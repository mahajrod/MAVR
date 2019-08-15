#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import FileRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file")
parser.add_argument("-c", "--column_number", action="store", dest="column_number", type=int, required=True,
                    help="Column number with values to be replaced(0-based)")
parser.add_argument("-m", "--comment_prefix", action="store", dest="comment_prefix", default="#",
                    help="Prefix of comment lines. Default: #")
parser.add_argument("-s", "--synonym_file", action="store", dest="syn_file", required=True,
                    help="Synonym file")
parser.add_argument("-e", "--separator", action="store", dest="separator", default="\t",
                    help="Separator used in input file")
parser.add_argument("-o", "--output", action="store", dest="output", default=None,
                    help="Output file")

args = parser.parse_args()

FileRoutines.replace_column_value_by_syn(args.input, args.syn_file, args.output, column=args.column_number,
                                         comment_prefix=args.comment_prefix, separator=args.separator,
                                         syn_header=False, syn_separator="\t",
                                         syn_key_index=0, syn_value_index=1, syn_comment_prefix=None)
