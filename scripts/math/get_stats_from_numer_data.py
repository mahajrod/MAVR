#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

from Routines import MathRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with data")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Prefix of output file")
parser.add_argument("-n", "--min_value", action="store", dest="min_length", type=float,
                    help="Minimum value to use. Default - not set")
parser.add_argument("-x", "--max_value", action="store", dest="max_length", type=float,
                    help="Maximum value to use. Default - not set")
parser.add_argument("-c", "--comments_prefix", action="store", dest="comments_prefix", default="#",
                    help="Comments prefix. Default - '#'")
parser.add_argument("-d", "--delimiter", action="store", dest="delimiter",
                    help="Value delimiter. Default - any space symbol")

args = parser.parse_args()
MathRoutines.get_stats_from_file(args.input, minimum=args.min_value, maximum=args.max_value, dtype=float,
                                 comments=args.comments_prefix, delimiter=args.delimiter, converters=None, skiprows=0,
                                 usecols=None, unpack=False, ndmin=0, output_file=args.output, verbose=True)
