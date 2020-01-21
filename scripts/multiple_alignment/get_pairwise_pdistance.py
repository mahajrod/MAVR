#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.Parsers.Sequence import CollectionSequence
from RouToolPa.Routines import MultipleAlignmentRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with alignment or sequences of equal length")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="File to write table of p-distances. Default: stdout")

args = parser.parse_args()

collection_sequence = CollectionSequence(in_file=args.input, parsing_mode="parse")

MultipleAlignmentRoutines.get_pairwise_pdist(collection_sequence.records, outfile=args.output)
