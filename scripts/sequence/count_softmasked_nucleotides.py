#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import SequenceRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input file with softmasked sequences")
parser.add_argument("-o", "--output_file_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of input and output files. Allowed formats genbank, fasta(default)")
parser.add_argument("-p", "--print_stats", action="store_true", dest="print_stats",
                    help="Print stats. Default: False")

args = parser.parse_args()

SequenceRoutines.count_softmasked_nucleotides_from_file(args.input, args.output_prefix, verbose=args.print_stats, parsing_mode="parse",
                                                        format=args.format, index_file=None)
