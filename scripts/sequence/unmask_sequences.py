#!/usr/bin/env python
__author__ = 'mahajrod'
import os
import sys
import argparse
from RouToolPa.Parsers.Sequence import CollectionSequence


parser = argparse.ArgumentParser()

parser.add_argument("-f", "--format", action="store", dest="format",
                    default="fasta", help="format of file with sequences - default: fasta.")
parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with sequences")
parser.add_argument("-o", "--output", action="store", dest="output", default="out.t",
                    help="Output file - default: out.t.")
parser.add_argument("-m", "--mode", action="store", dest="mode", default="parse",
                    help="Parsing mode. Allowed: parse(default), generator")

args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

seq_collection = CollectionSequence(in_file=args.input, format=args.format,
                                    parsing_mode=args.mode)
seq_collection.unmask(in_place=True)
seq_collection.write(args.output)
