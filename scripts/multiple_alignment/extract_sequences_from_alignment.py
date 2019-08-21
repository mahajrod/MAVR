#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import MultipleAlignmentRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with alignment")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of the input alignment. Default: fasta")
parser.add_argument("-a", "--output_format", action="store", dest="output_format", default="fasta",
                    help="Format of the output sequences. Default: fasta")
parser.add_argument("-g", "--gap_symbol", action="store", dest="gap_symbol", default="-",
                    help="Gap symbol. Default: -")

args = parser.parse_args()

MultipleAlignmentRoutines.extract_sequences_from_alignment(args.input,
                                                           args.output,
                                                           alignment_format=args.format,
                                                           output_format=args.output_format,
                                                           gap_symbol=args.gap_symbol)
