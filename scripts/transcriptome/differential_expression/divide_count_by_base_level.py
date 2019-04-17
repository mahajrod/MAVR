#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import ExpressionRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with counts")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file with counts")
parser.add_argument("-b", "--base_level", action="store", dest="base_level", required=True,
                    help="Base level id (must be present in header)")
parser.add_argument("-s", "--secondary_base_level", action="store", dest="secondary_base_level",
                    help="Secondary base level id (must be present in header if set)")

args = parser.parse_args()

ExpressionRoutines.divide_counts_by_base_level(args.input, args.output, args.base_level,
                                               separator="\t",
                                               secondary_base_lvl=args.secondary_base_level)
