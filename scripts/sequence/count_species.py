#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import SequenceRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input genbank file")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file with species counts")
parser.add_argument("-f", "--format", action="store", dest="format", default="genbank",
                    help="Format of input file. Default - genbank ")

args = parser.parse_args()

SequenceRoutines.count_species_from_file(args.input, format=args.format,
                                         output_filename="count_species.count")
