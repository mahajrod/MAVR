#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Tools.Bedtools import GenomeCov

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input",
                    help="Input collapsed coverage file")
parser.add_argument("-o", "--output", action="store", dest="output",
                    help="Output file with stats")

args = parser.parse_args()

GenomeCov.analyze_collapsed_coverage_file(args.input, args.output)
