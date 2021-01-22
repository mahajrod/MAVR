#!/usr/bin/env python
__author__ = 'mahajrod'

import argparse
from RouToolPa.Tools.Alignment import Mosdepth

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Per base bed coverage file")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-w", "--window_size", action="store", dest="window_size", default=100000, type=int,
                    help="Size of windows. Default: 100 000")
parser.add_argument("-s", "--window_step", action="store", dest="window_step", default=None, type=int,
                    help="Step of windows. Default: window size")
parser.add_argument("-b", "--buffer_size", action="store", dest="buffer_size", default=10000000, type=int,
                    help="Buffer size in bytes for reading file. Default: 10 000 000")
args = parser.parse_args()

Mosdepth.get_coverage_stats_in_windows(args.input, args.window_size, args.output_prefix, window_step=args.window_step,
                                       buffering=args.buffer_size)
