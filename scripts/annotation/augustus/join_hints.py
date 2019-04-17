#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Annotation import AUGUSTUS

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Comma-separated list of input files with hints in gff")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file with merged hints")


args = parser.parse_args()

AUGUSTUS.join_multiple_hints(args.input, args.output)
