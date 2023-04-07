#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse

import pandas as pd

from RouToolPa.Parsers.BED import CollectionBED

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input bed file. Default: stdin")
parser.add_argument("-r", "--right_flank", action="store", dest="right_flank", default=0, type=int,
                    help="Length of right flank to add. Default: 0")
parser.add_argument("-l", "--left_flank", action="store", dest="left_flank", default=0, type=int,
                    help="Length of left flank to add. Zero border is controlled. Default: 0")
parser.add_argument("-e", "--length_file", action="store", dest="length_file", default=None,
                    help="File with lengths of scaffolds. Default: not set")
parser.add_argument("-s", "--scaffold_column", action="store", dest="scaffold_column", default=0,  type=int,
                    help="Number (0-based) of column with scaffold ids. Default: 0")
parser.add_argument("-c", "--length_column", action="store", dest="length_column", default=1,  type=int,
                    help="Number (0-based) of column with lengths of scaffolds. Default: 1")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output bed file with flanks. Default: stdout")

args = parser.parse_args()
length_df = pd.read_csv(args.length_file, sep='\t', header=None, index_col=args.scaffold_column,
                        usecols=(args.scaffold_column, args.length_column)) if args.length_file else None
if length_df:
    length_df.columns = pd.Index(["length", ])
bed_col = CollectionBED(in_file=args.input, parsing_mode="all")

bed_col.add_flanks(left_flank=args.left_flank, right_flank=args.right_flank, length_df=None,
                   length_df_column="length", inplace=True)
bed_col.write(args.output)

