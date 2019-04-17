#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import MathRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with data table")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file")
parser.add_argument("-c", "--comments_prefix", action="store", dest="comments_prefix", default="#",
                    help="Comments prefix. Default - '#'")
parser.add_argument("-s", "--separator", action="store", dest="separator", default='\t',
                    help="Separator between columns in table. Default: '\\t' ")

args = parser.parse_args()
MathRoutines.get_per_row_stats_for_table_from_file(args.input, args.output, rownames=True, column_names=True,
                                                   separator=args.separator, comments_prefix=args.comments_prefix)
