#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.Routines import FileRoutines



parser = argparse.ArgumentParser()

parser.add_argument("-f", "--file_list", action="store", dest="file_list", required=True,
                    type=FileRoutines.make_list_of_path_to_files_from_string,
                    help="Comma-separated list of files or directorie with files to combine")

parser.add_argument("-c", "--column_index_list", action="store", dest="column_index_list", required=True,
                    type=lambda s: map(lambda r: map(int, r.split(",")), s.split(":")),
                    help="Colon-separated list of comma-separated list of  0-based column indexes to use for "
                         "file merging. Example: 0,1;0,1;0,1 if you want to merge three files by 0 and 1st columns")
parser.add_argument("-e", "--header", action="store_true", dest="header",
                    help="Files contain headers. Default - False")
parser.add_argument("-n", "--column_names_list", action="store", dest="column_names_list",
                    type=lambda s: s.split(","),
                    help="Comma-separated list of new names for columns used for file merging."
                         "Necessary if -h/--header option is specified, otherwise it is ignored")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output file with allowed ids. Default - stdout")

args = parser.parse_args()


FileRoutines.merge_files_by_columns(args.file_list, args.column_index_list, args.output, separator="\t",
                                    column_names_list=args.column_names_list, comment_prefix="#", header=args.header)
