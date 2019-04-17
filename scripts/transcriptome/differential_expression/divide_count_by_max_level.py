#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import ExpressionRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with counts")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file with counts")

args = parser.parse_args()

ExpressionRoutines.divide_counts_by_max_level(args.input, args.output,
                                              separator="\t",)
