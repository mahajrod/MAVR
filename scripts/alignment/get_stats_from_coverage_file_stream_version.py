#!/usr/bin/env python
__author__ = 'mahajrod'

import argparse
from RouToolPa.Tools.Bedtools import GenomeCov

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Coverage file")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file with stats")
parser.add_argument("-s", "--scaffold_column", action="store", dest="scaffold_column", default=0, type=int,
                    help="Column(0-based) with scaffold ids. Default: 0")
parser.add_argument("-c", "--coverage_column", action="store", dest="coverage_column", default=1, type=int,
                    help="Column(0-based) with per base coverage values. Default: 1")
args = parser.parse_args()

GenomeCov.get_stats_from_coverage_file_stream_version(args.input, args.output, verbose=True,
                                                      scaffold_column=args.scaffold_column,
                                                      coverage_column=args.coverage_column,
                                                      separator="\t")
