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
                    help="Parsing mode. Allowed: parse(default), index, index_db")

args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

coll_seq = CollectionSequence(in_file=args.input, parsing_mode="generator", get_stats=True)
coll_seq.seq_lengths.to_csv(args.output, header=False)

#record_dict = SequenceRoutines.parse_seq_file(args.input, args.mode, format=args.format, index_file="temp_index.idx")
#lengths_dict = SequenceRoutines.get_lengths(record_dict, out_file=out_fd)

print("Longest sequence: %i" % max(coll_seq.seq_lengths["length"]))
print("Shortest sequence: %i" % min(coll_seq.seq_lengths["length"]))
print("Total length: %i" % sum(coll_seq.seq_lengths["length"]))

#if args.mode == "index_db":
#    os.remove("temp_index.idx")
