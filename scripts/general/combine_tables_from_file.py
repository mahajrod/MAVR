#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import FileRoutines
from RouToolPa.Collections.General import TwoLvlDict



parser = argparse.ArgumentParser()

parser.add_argument("-f", "--files", action="store", dest="files", required=True,
                    type=FileRoutines.make_list_of_path_to_files_from_string,
                    help="Comma-separated list of files/directories with tables")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file with combined table.")
parser.add_argument("-a", "--absent_symbol", action="store", dest="absent_symbol", default=".",
                    help="Symbol to be treated as absent value")
parser.add_argument("-v", "--split_values", action="store_true", dest="split_values",
                    help="Split values. Default: False")
parser.add_argument("-s", "--value_separator", action="store", dest="value_separator", default=",'",
                    help="Value separator. Default: ','")
parser.add_argument("-g", "--ignore_value_repeats", action="store_true", dest="ignore_value_repeats",
                    help="Ignore repeats of values(i.e values that corresponds to same fl_key and sl_key) "
                         "and don't raise exception. If yes value from first entry is stored. Default: False")

args = parser.parse_args()

combined_table = TwoLvlDict(input_file=args.files, absent_symbol=args.absent_symbol,
                            split_values=args.split_values, value_sep=args.value_separator,
                            ignore_value_repeats=args.ignore_value_repeats)
#print combined_table
combined_table.write(args.output, absent_symbol=args.absent_symbol, close_after_if_file_object=False, sort=False)
