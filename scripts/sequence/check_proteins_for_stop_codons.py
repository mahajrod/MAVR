#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import SequenceRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input file with protein sequences")
parser.add_argument("-o", "--output_file_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of input and output files. Allowed formats genbank, fasta(default)")
parser.add_argument("-s", "--stop_codons", action="store", dest="stop_codons",
                    default=("*", "."),
                    type=lambda s: set(s.split(",")),
                    help="Comma-separated list of stop codon symbols. Default - '.', '*'")
args = parser.parse_args()

SequenceRoutines.check_proteins_for_stop_codons_from_file(args.input, args.output_prefix,
                                                          stop_codon_symbol_set=args.stop_codons,
                                                          format=args.format)
