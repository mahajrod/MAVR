#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.Routines import FileRoutines



parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    type=FileRoutines.make_list_of_path_to_files_from_string,
                    help="Comma-separated list of files/directories")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output merged file")
parser.add_argument("-e", "--header_prefix", action="store", dest="header_prefix", default="#",
                    help="Header prefix")
parser.add_argument("-s", "--sorting_options", action="store", dest="sorting_options",
                    help="Sorting options for sort utility. Default: not set, i.e. no sort")


args = parser.parse_args()

FileRoutines.combine_files_with_header(args.input,
                                       args.output,
                                       header_prefix=args.header_prefix,
                                       sorting_options=args.sorting_options)
