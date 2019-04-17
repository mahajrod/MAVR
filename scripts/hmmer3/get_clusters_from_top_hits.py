#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.HMMER import HMMER3

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input",
                    help="File with HMMER top hits")
parser.add_argument("-o", "--output", action="store", dest="output",
                    help="Output file")
args = parser.parse_args()

HMMER3.get_families_from_top_hits(args.input, args.output)
