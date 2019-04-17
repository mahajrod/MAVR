#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Tools.Annotation import AUGUSTUS


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input evidence file")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="File to write filtered evidence file")
parser.add_argument("-m", "--min_fraction", action="store", dest="min_fraction", default=0, type=float,
                    help="Minimum fraction of transcript supported by hints")

args = parser.parse_args()

AUGUSTUS.extract_longest_isoforms(args.input, args.output,
                                  minimum_supported_fraction=args.min_fraction)
