#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import sys
import argparse
from RouToolPa.Collections.General import SynDict
from RouToolPa.Routines.File import make_list_of_path_to_files



def make_list_of_path_to_files_from_comma_sep_string(string):
    return make_list_of_path_to_files(string.split(","))

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    type=make_list_of_path_to_files_from_comma_sep_string,
                    help="Comma-separated list of fam files or directories with them")

parser.add_argument("-k", "--key_column_index", action="store", dest="key_column_index", type=int, default=0,
                    help="Index of key column in synonym file. Default: 0")
parser.add_argument("-v", "--value_column_index", action="store", dest="value_column_index", type=int, default=1,
                    help="Index of value column in synonym file.Default: 1")
parser.add_argument("-a", "--allow_repeats_of_key", action="store_true", dest="allow_repeats_of_key", default=False,
                    help="Allow repeats of keys. Default: False")

parser.add_argument("-o", "--output_file", action="store", dest="output", default="stdout",
                    help="Output file")


args = parser.parse_args()


out_fd = sys.stdout if args.input == "stdout" else open(args.output, "w")
family_dict = SynDict()

for filename in args.input:
    fam_dict = SynDict()
    fam_dict.read(filename, split_values=True, allow_repeats_of_key=args.allow_repeats_of_key,
                  key_index=args.key_column_index, value_index=args.value_column_index)
    for family in fam_dict:
        if family not in family_dict:
            family_dict[family] = fam_dict[family]
        else:
            family_dict[family] += fam_dict[family]

family_dict.write(args.output, splited_values=True)

if args.output != "stdout":
    out_fd.close()
