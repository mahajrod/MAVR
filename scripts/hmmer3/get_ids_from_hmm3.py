#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.HMMER import HMMER3

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_hmm", action="store", dest="input",
                    help="Input hmm3 file")
parser.add_argument("-o", "--output_file", action="store", dest="output", default="stdout",
                    help="Output file with ids. Default: stdout")

args = parser.parse_args()

HMMER3.get_ids_from_hmm3(args.input, ids_file=None if args.output == 'stdout' else args.output, return_ids_list=False)
