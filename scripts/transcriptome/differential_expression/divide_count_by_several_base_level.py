#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import ExpressionRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with counts")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-b", "--base_levels", action="store", dest="base_levels", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of base levels")
parser.add_argument("-m", "--max_ratio_to_base_level", action="store", dest="max_ratio_to_base_level",
                    type=float, default=0.5,
                    help="Maximum ratio between nonbase level and max(base level) to retain in ratio file. "
                         "Default: 0.5")

args = parser.parse_args()

ExpressionRoutines.divide_counts_by_several_base_level(args.input, args.output_prefix, args.base_levels,
                                                       separator="\t", verbose=True,
                                                       max_ratio_to_base_lvl=args.max_ratio_to_base_level)
