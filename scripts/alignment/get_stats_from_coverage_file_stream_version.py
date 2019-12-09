#!/usr/bin/env python
__author__ = 'mahajrod'

import argparse
from RouToolPa.Tools.Bedtools import GenomeCov

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Coverage file")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file with stats")

args = parser.parse_args()

GenomeCov.get_stats_from_coverage_file_stream_version(args.input, args.output, verbose=True,
                                                      scaffold_column=0, coverage_column=1,
                                                      separator="\t")
