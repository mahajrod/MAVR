#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.Parsers.Sequence import CollectionSequence
from RouToolPa.Routines import MultipleAlignmentRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with alignment or sequences of equal length")
parser.add_argument("-p", "--output_prefix", action="store", dest="output_prefix",
                    help="Prefix of output_files")
parser.add_argument("-m", "--max_distance", action="store", dest="min_distance", type=int,
                    help="Maximal distance between sequences to treat them as close. Default: not set")
parser.add_argument("-a", "--any_symbol", action="store", dest="any_symbol", default=None,
                    help="Symbol to treated as equal to any other symbol. Use N for nucleotides and X - for proteins. "
                         "Default: not set")

args = parser.parse_args()

collection_sequence = CollectionSequence(in_file=args.input, parsing_mode="parse")

MultipleAlignmentRoutines.get_pairwise_hamming(collection_sequence.records, output_prefix=args.output_prefix,
                                               min_distance=args.min_distance,
                                               any_symbol=args.any_symbol)
