#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import SequenceRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input file with cds sequences")
parser.add_argument("-o", "--output_file_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of input and output files. Allowed formats genbank, fasta(default)")
parser.add_argument("-g", "--genetic_code", action="store", dest="genetic_code",
                    default=1, type=int,
                    help="Number of NCBI genetic code table to use. Default - 1(universal) ")
args = parser.parse_args()

SequenceRoutines.check_cds_for_stop_codons_from_file(args.input, args.output_prefix,
                                                     genetic_code_table=args.genetic_code, format=args.format)
