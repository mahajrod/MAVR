#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.Routines import AlignmentRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-m", "--mask_file", action="store", dest="mask_file", required=True,
                    help="File with per base mask")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file with collapsed(from bases to regions) mask")

parser.add_argument("--scaf_col", action="store", dest="scaf_col", type=int,
                    default=0,
                    help="Index of scaffold column in coverage file (0-based). Default: 0")
parser.add_argument("--pos_col", action="store", dest="pos_col", type=int,
                    default=1,
                    help="Index of position column in coverage file (0-based). Default: 1")
parser.add_argument("-z", "--0-based", action="store_true", dest="zero_based_output",
                    help="Write output in 0-based format and python notation, "
                         "i.e last position is not included in the interval."
                         "Required for bed format. Default: False")
args = parser.parse_args()

AlignmentRoutines.collapse_per_base_coverage_mask(args.mask_file, args.output,
                                                  scaffold_column=args.scaf_col,
                                                  position_column=args.pos_col,
                                                  zero_based_output=args.zero_based_output,
                                                  comments_prefix="#")
