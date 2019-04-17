#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.Routines import FileRoutines



parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output file")
parser.add_argument("-k", "--key_column", action="store", dest="key_column", type=int, required=True,
                    help="Column to be used as key(0-based). ")
parser.add_argument("-y", "--syn_file", action="store", dest="syn_file", required=True,
                    help="Synonym file with new column values")
parser.add_argument("-s", "--separator", action="store", dest="separator", default="\t",
                    help="Separator used in input file. Default: TAB")

parser.add_argument("-n", "--new_column_name", action="store", dest="new_column_name",
                    required=True,
                    help="Name of new column")

args = parser.parse_args()

FileRoutines.add_add_new_column_by_key_column(args.input,
                                              args.syn_file,
                                              args.key_column,
                                              args.output,
                                              new_column_name=args.new_column_name,
                                              separator=args.separator,
                                              absent_value=".")
