#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import SequenceRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Comma separated list of genbank files/directories")
parser.add_argument("-o", "--output_file_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-f", "--format", action="store", dest="format", default="genbank",
                    help="Format of input and output file. Allowed formats genbank(default), fasta")

args = parser.parse_args()

SequenceRoutines.get_random_species_genomes_from_genbank_file(args.input, args.output_prefix,
                                                              output_type=args.format)
