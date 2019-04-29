#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.WGA import LAST


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_tab_file", action="store", dest="input_tab_file", required=True,
                    help="Input TAB file")
parser.add_argument("-o", "--output_tab_file", action="store", dest="output_tab_file", required=True,
                    help="Output TAB file")
parser.add_argument("-l", "--min_len", action="store", dest="min_len", required=True, type=int,
                    help="Minimum length of hit to retain.")
parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                    help="Verbose output")
parser.add_argument("-m", "--mode", action="store", dest="mode", default="target",
                    help="Filtering mode. Allowed: target (default), query")

args = parser.parse_args()


LAST.filter_tab_output_by_length(args.input_tab_file, args.min_len, args.output_tab_file,
                                 mode=args.mode, verbose=args.verbose)
