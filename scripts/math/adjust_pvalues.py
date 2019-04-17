#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import MathRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file p-values to be adjusted")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file to write results")
parser.add_argument("-n", "--pvalue_column", action="store", dest="pvalue_column", type=int, required=True,
                    help="0-based number of column with raw p-values")
parser.add_argument("-m", "--max", action="store", dest="max", type=float,
                    help="Maximum value to use. Default - not set")
parser.add_argument("-c", "--comments_prefix", action="store", dest="comments_prefix", default="#",
                    help="Comments prefix. Default - '#'")
parser.add_argument("-s", "--separator", action="store", dest="separator", default="\t",
                    help="Separator used in input file. Default - \\t")

args = parser.parse_args()
MathRoutines.adjust_pvalues_from_file(args.input, args.pvalue_column, args.output,
                                      header=True, comments_prefix=args.comments_prefix, separator=args.separator)
