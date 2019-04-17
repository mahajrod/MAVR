#!/usr/bin/env python
__author__ = 'mahajrod'
import os
import sys
import argparse
from RouToolPa.Routines.Sequence import SequenceRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-f", "--format", action="store", dest="format",
                    default="fasta", help="format of file with sequences - default: fasta.")
parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with sequences")
parser.add_argument("-o", "--output", action="store", dest="output", default="out.t",
                    help="Output file - default: out.t.")
parser.add_argument("-m", "--mode", action="store", dest="mode", default="parse",
                    help="Parsing mode. Allowed: parse(default), index, index_db")

args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

record_dict = SequenceRoutines.parse_seq_file(args.input, args.mode, format=args.format, index_file="temp_index.idx")
lengths_dict = SequenceRoutines.get_lengths(record_dict, out_file=out_fd)

print("Longest sequence: %i" % max(lengths_dict.values()))
print("Shortest sequence: %i" % min(lengths_dict.values()))
print("Total length: %i" % sum(lengths_dict.values()))

if args.mode == "index_db":
    os.remove("temp_index.idx")
