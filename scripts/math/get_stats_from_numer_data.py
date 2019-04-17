#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import MathRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with data")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file")
parser.add_argument("-n", "--min", action="store", dest="min", type=float,
                    help="Minimum value to use. Default - not set")
parser.add_argument("-m", "--max", action="store", dest="max", type=float,
                    help="Maximum value to use. Default - not set")
parser.add_argument("-c", "--comments_prefix", action="store", dest="comments_prefix", default="#",
                    help="Comments prefix. Default - '#'")
parser.add_argument("-d", "--delimiter", action="store", dest="delimiter",
                    help="Value delimiter. Default - any space symbol")
parser.add_argument("-l", "--columns", action="store", dest="columns", type=lambda s: map(int, s.split(",")),
                    help="Comma-separated list of columns(0-based) with data. Default - all")

args = parser.parse_args()

MathRoutines.get_stats_from_file(args.input, minimum=args.min, maximum=args.max, dtype=float,
                                 comments=args.comments_prefix, delimiter=args.delimiter, converters=None, skiprows=0,
                                 usecols=args.columns, unpack=False, ndmin=0, output_file=args.output, verbose=True)
