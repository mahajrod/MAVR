#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from Routines.File import tsv_extract_by_column_value


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input_file",
                    help="Input file")
parser.add_argument("-c", "--column_number", action="store", dest="column_number", type=int,
                    help="Column number")
parser.add_argument("-v", "--values", action="store", dest="values",
                    help="Lines with this values in corresponding column will be extracted. Values should be separated by commas")
parser.add_argument("-s", "--separator", action="store", dest="separator", default="\t",
                    help="Separator used in input file")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", default=None,
                    help="Column number")
parser.add_argument("-e", "--header", action="store_true", dest="header", default=False,
                    help="Set if header is present in input file")

args = parser.parse_args()

args.values = args.values.split(",")
tsv_extract_by_column_value(args.input_file, args.column_number, args.values, separator=args.separator,
                             header=args.header, outfile_prefix=args.output_prefix)