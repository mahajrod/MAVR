#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.Collections.General import SynDict



parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input file")
parser.add_argument("-o", "--output_file", action="store", dest="output", default="stdout",
                    help="Output file with collapsed strings")
parser.add_argument("-c", "--column_separator", action="store", dest="column_separator", default="\t",
                    help="Column separator. Default: '\\t'")
parser.add_argument("-v", "--value_separator", action="store", dest="value_separator", default=",",
                    help="Value separator. Default: ','")
parser.add_argument("-k", "--key_column", action="store", dest="key_column", default=0, type=int,
                    help="Column to be used as key(0-based). Default: 0")
parser.add_argument("-a", "--value_column", action="store", dest="value_column", default=1, type=int,
                    help="Column to be used as value(0-based). Default: 1")
parser.add_argument("-m", "--comments_prefix", action="store", dest="comments_prefix", default="#",
                    help="Prefix of strings(comments) to be ignored. Default: #")
parser.add_argument("-r", "--remove_value_repeats", action="store_true", dest="remove_value_repeats",
                    help="Remove repeats of values")
args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

syn_dict = SynDict()
syn_dict.read(args.input, header=False, separator=args.column_separator, allow_repeats_of_key=True,
              split_values=True, values_separator=args.value_separator,
              key_index=args.key_column, value_index=args.value_column,
              comments_prefix=args.comments_prefix)

if args.remove_value_repeats:
    collapsed_dict = syn_dict.remove_value_repeats()
    collapsed_dict.write(out_fd, splited_values=True, values_separator=args.value_separator,
                         close_after_if_file_object=True)
else:
    syn_dict.write(out_fd, splited_values=True, values_separator=args.value_separator,
                   close_after_if_file_object=True)
#out_fd.close()