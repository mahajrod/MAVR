#!/usr/bin/env python
__author__ = 'mahajrod'
import sys
import argparse
from RouToolPa.Parsers.Sequence import CollectionSequence

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Fasta file with sequences. Default: stdin")

parser.add_argument("-o", "--output_directory", action="store", dest="out_dir",
                    help="Directory to write output files")


args = parser.parse_args()

sequence_collection = CollectionSequence(in_file=args.input, parsing_mode='parse')
sequence_collection.write_splited(args.out_dir)
