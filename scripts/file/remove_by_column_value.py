#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Collections.General import IdList
from RouToolPa.Routines import FileRoutines




parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input_file",
                    help="Input file")
parser.add_argument("-c", "--column_number", action="store", dest="column_number", type=int,
                    help="Column number(0-based)")
parser.add_argument("-v", "--values", action="store", dest="values",
                    help="Lines with this values in corresponding column will be removed."
                         "Values should be separated by commas")
parser.add_argument("-f", "--values_file", action="store", dest="values_file",
                    help="File with values. Ignored if -v/--values is set")
parser.add_argument("-s", "--separator", action="store", dest="separator", default="\t",
                    help="Separator used in input file. Default: '\t'")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", default=None,
                    help="Output_prefix")
parser.add_argument("-e", "--header", action="store_true", dest="header", default=False,
                    help="Set if header is present in input file")

args = parser.parse_args()

if (not args.values) and (not args.values_file):
    raise ValueError("Neither -v/--values nor -f/--values_file was set")

values = args.values.split(",") if args.values else IdList(filename=args.values_file)

FileRoutines.tsv_remove_by_column_value(args.input_file, args.column_number, values, separator=args.separator,
                                        header=args.header, outfile_prefix=args.output_prefix)
