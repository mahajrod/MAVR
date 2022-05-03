#!/usr/bin/env python
__author__ = 'mahajrod'
import os
import sys
import argparse
from RouToolPa.Parsers.Sequence import CollectionSequence
from RouToolPa.GeneralRoutines import FileRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-f", "--format", action="store", dest="format",
                    default="fasta", help="format of file with sequences - default: fasta.")
parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with sequences")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    type=lambda s: FileRoutines.metaopen(s, "w"),
                    help="Output file - default: stdout")
parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                    help="Write some stats to stderr. Default: False")
args = parser.parse_args()

coll_seq = CollectionSequence(in_file=args.input, parsing_mode="generator", get_stats=True)
coll_seq.seq_lengths.to_csv(args.output, header=False, sep="\t")

if args.verbose:
    sys.stderr.write("Longest sequence: %i\n" % max(coll_seq.seq_lengths["length"]))
    sys.stderr.write("Shortest sequence: %i\n" % min(coll_seq.seq_lengths["length"]))
    sys.stderr.write("Total length: %i\n" % sum(coll_seq.seq_lengths["length"]))
