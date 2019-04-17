#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.Routines import FileRoutines



parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input",  default=sys.stdin,
                    help="Input file")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output file")
parser.add_argument("-k", "--column_number", action="store", dest="column_number", type=int, required=True,
                    help="Column to be get values from (0-based).")
parser.add_argument("-s", "--separator", action="store", dest="separator", default="\t",
                    help="Separator used in input file. Default: TAB")

args = parser.parse_args()

FileRoutines.get_column_value_set_from_file(args.input, args.column_number, output_file=args.output,
                                            separator=args.separator, comments_prefix="#", verbose=False)
