#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import MultipleAlignmentRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with codon alignment")

parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of thr input alignment. Default: fasta")

args = parser.parse_args()

MultipleAlignmentRoutines.extract_codon_positions_from_file(args.input, args.output_prefix, format=args.format)
