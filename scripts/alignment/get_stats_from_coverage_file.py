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

GenomeCov.get_coverage_stats(args.input, args.output, verbose=True)
