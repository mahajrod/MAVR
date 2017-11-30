#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from BConverters import SequenceConverters

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input genbank file")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output fasta file")
parser.add_argument("-p", "--parsing_mode", action="store", dest="parsing_mode",
                    default="parse", help="Parsing mode. Allowed: parse(default), index, index_db")

args = parser.parse_args()

SequenceConverters.convert_sequences(args.input, "genbank", args.output, "fasta",
                                     parsing_mode=args.parsing_mode, input_index="tmp.idx")

