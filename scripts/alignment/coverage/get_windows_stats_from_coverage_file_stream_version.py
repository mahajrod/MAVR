#!/usr/bin/env python
__author__ = 'mahajrod'

import argparse
from RouToolPa.Parsers.AGP import CollectionAGP

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Coverage file")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file")
parser.add_argument("-w", "--window_size", action="store", dest="window_size", default=100000, type=int,
                    help="Size of windows. Default: 100 000")
parser.add_argument("-s", "--window_step", action="store", dest="window_step", default=None, type=int,
                    help="Step of windows. Default: window size")
parser.add_argument("-b", "--buffer_size", action="store", dest="buffer_size", default=10000000, type=int,
                    help="Buffer size in bytes for reading file. Default: 10 000 000")
args = parser.parse_args()

GenomeCov.get_coverage_stats_in_windows(args.input, args.window_size, args.output, window_step=args.window_step,
                                        buffering=args.buffer_size)
