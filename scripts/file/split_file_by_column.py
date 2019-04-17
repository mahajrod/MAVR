#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import FileRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input_file", required=True,
                    help="Input file")
parser.add_argument("-c", "--column_number", action="store", dest="column_number", type=int, required=True,
                    help="Column number to use for split(0-based)")
parser.add_argument("-s", "--separator", action="store", dest="separator", default="\t",
                    help="Separator used in input file")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", default=None,
                    help="Prefix of output files")
parser.add_argument("-e", "--header", action="store_true", dest="header", default=False,
                    help="Set if header is present in input file")
parser.add_argument("-u", "--use_column_value_as_prefix", action="store_true", dest="use_column_value_as_prefix",
                    default=False,
                    help="Use column value as prefix for output files")
parser.add_argument("-r", "--sorted_input", action="store_true", dest="sorted_input",
                    default=False,
                    help="Input file is sorted. Do it to reduce number of simultaneously opened files")


args = parser.parse_args()

FileRoutines.split_by_column(args.input_file, args.column_number, separator=args.separator,
                             header=args.header, outfile_prefix=args.output_prefix,
                             use_column_value_as_prefix=args.use_column_value_as_prefix, sorted_input=args.sorted_input)
