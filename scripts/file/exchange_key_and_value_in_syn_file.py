#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.Collections.General import SynDict




parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input syn file")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output syn file")
parser.add_argument("-s", "--separator", action="store", dest="separator", default="\t",
                    help="Separator used in input file. Default: TAB")
parser.add_argument("-l", "--split_values", action="store_true", dest="split_values", default=False,
                    help="Split values in input file")
parser.add_argument("-e", "--value_separator", action="store", dest="value_separator", default=",",
                    help="Value separator in input file. Default: ,")

args = parser.parse_args()

input_dict = SynDict(filename=args.input, separator=args.separator, split_values=args.split_values,
                     values_separator=args.value_separator)

output_dict = input_dict.exchange_key_and_value()
output_dict.write(args.output, splited_values=True)
