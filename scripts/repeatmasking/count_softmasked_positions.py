#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.Parsers.Sequence import CollectionSequence

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_fasta", action="store", dest="input_fasta", default=sys.stdin,
                    help="Input fasta file with sequences. Default: stdin")
parser.add_argument("-o", "--output", action="store", dest="output", default=None,
                    help="Output file with softmasked position counts per input sequence."
                         " If output file is not stdout, a total count will be printed. "
                         " Default: stdout")
args = parser.parse_args()


CollectionSequence(in_file=args.input_fasta).count_masked_positions(out_file=sys.stdout if args.output is None else args.output,
                                                                    verbose=True if args.output is not None else False)
