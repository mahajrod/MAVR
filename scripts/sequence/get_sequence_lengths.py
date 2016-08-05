#!/usr/bin/env python
__author__ = 'mahajrod'
import os
import sys
import argparse

from Bio import SeqIO

from Routines.Sequence import SequenceRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-f", "--format", action="store", dest="format",
                    default="fasta", help="format of file with sequences - default: fasta.")
parser.add_argument("-i", "--input", action="store", dest="input",
                    help="file with sequences")
parser.add_argument("-o", "--output", action="store", dest="output",
                    default="out.t", help="output file - default: out.t.")

args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

record_dict = SeqIO.index_db("temp_index.idx", [args.input], format=args.format)
lengths_dict = SequenceRoutines.get_lengths(record_dict, out_file=out_fd)
print("Longest sequence: %i" % max(lengths_dict.values()))
print("Shortest sequence: %i" % min(lengths_dict.values()))
print("Total length: %i" % sum(lengths_dict.values()))
os.remove("temp_index.idx")