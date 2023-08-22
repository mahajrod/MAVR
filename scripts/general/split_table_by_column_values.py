#!/usr/bin/env python
__author__ = 'Sergei Kliver'
from pathlib import Path
from io import TextIOWrapper
import sys
import argparse
import pandas as pd
from RouToolPa.GeneralRoutines import FileRoutines
from RouToolPa.Collections.General import TwoLvlDict


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input file with table.  Default: stdin")
parser.add_argument("-c", "--column", action="store", dest="column", required=True,
                    help="Column number (zero-based) or column name if header is present. "
                         "Script will first check if value is convertable to int, and if yes will use it as "
                         "column index. If not script will use value as column name. "
                         "In such case if  --header is not set an exception will be risen. Required.")
parser.add_argument("-k", "--keep_column", action="store_true", dest="keep_column", default=False,
                    help="Keep column used for split in the output. Default: False")
parser.add_argument("-a", "--header", action="store_true", dest="header", default=None,
                    help="Header is present in input. Default: False")
parser.add_argument("-s", "--input_separator", action="store", dest="input_separator", default="\t",
                    help="Column separator in input file. Default: '\\t'")
parser.add_argument("-r", "--output_separator", action="store", dest="output_separator", default=None,
                    help="Column separator in output files. Default: same as --input_separator")
parser.add_argument("-e", "--output_extension", action="store", dest="output_extension", default=None,
                    help="Extension of output files (with leading dot). Required if reading input from stdin."
                         " If set when input is a file it will be used instead extension of the input files."
                         " Default: not set.")
parser.add_argument("-o", "--output_prefix", dest="output_prefix", required=True,
                    help="Output prefix. Extension for output files is detected automatically from the input file."
                         "Output files are names according the following template: "
                         " <output_prefix>.<column_value>.<extension>. Required.")

args = parser.parse_args()

output_separator = args.output_separator if args.output_separator is not None else args.input_separator
if isinstance(args.input, TextIOWrapper):
    if args.output_extension is None:
        raise ValueError("ERROR!!! Input is stdin and --output_extension is not set!")
    else:
        extension = args.output_extension
else:
    extension = Path(args.input).suffix if args.output_extension is None else args.output_extension

try:
    column_index = int(args.column)
except:
    column_index = None

if (column_index is None) and (args.header is None):
    raise ValueError("Value of --column is not convertable to int and --header is not set.")

data_df = pd.read_csv(args.input, sep=args.input_separator, header=0 if args.header is not None else None)
if args.keep_column:
    columns = data_df.columns
else:
    columns = list(data_df.columns)
    if column_index is None:
        columns.remove(args.column)
    else:
        columns.pop(column_index)
for value in (data_df[args.column] if column_index is None else data_df.iloc[:, column_index]).sort_values().unique():
    data_df.loc[(data_df[args.column] if column_index is None else data_df.iloc[:, column_index]) == value, columns].to_csv("{0}.{1}{2}".format(args.output_prefix,
                                                                                                                                   value,
                                                                                                                                   extension),
                                                                                                               sep=output_separator,
                                                                                                               header=True if args.header else False,
                                                                                                               index=False)
