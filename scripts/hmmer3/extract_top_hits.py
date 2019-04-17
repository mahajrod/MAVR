#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Tools.HMMER import HMMER3

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="File with hmm hits")
parser.add_argument("-p", "--output_prefix", action="store", dest="output_prefix",
                    help="Prefix of output file")
parser.add_argument("-a", "--parsing_mode", action="store", dest="parsing_mode", default="index_db",
                    help="Parsing mode for hmmer hits file. Allowed: parse, index, index_db(default)")

args = parser.parse_args()

HMMER3.extract_top_hits(args.input, args.output_prefix, parsing_mode=args.parsing_mode)
