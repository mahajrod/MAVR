#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.Routines import AlignmentRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-m", "--mask_file", action="store", dest="mask_file", required=True,
                    help="File with per base mask")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files with collapsed(from bases to regions) mask."
                         "Output is written in three formats: bed(0-based), tab, gff")
parser.add_argument("--scaf_col", action="store", dest="scaf_col", type=int,
                    default=0,
                    help="Index of scaffold column in coverage file (0-based). Default: 0")
parser.add_argument("--pos_col", action="store", dest="pos_col", type=int,
                    default=1,
                    help="Index of position column in coverage file (0-based). Default: 1")
parser.add_argument("--no_gzip", action="store_true", dest="no_gzip", default=False,
                    help="Don't gzip output files. Default: False")
args = parser.parse_args()

AlignmentRoutines.collapse_per_base_coverage_mask(args.mask_file, args.output_prefix,
                                                  scaffold_column=args.scaf_col,
                                                  position_column=args.pos_col,
                                                  comments_prefix="#",
                                                  gzip_output=not args.no_gzip)
