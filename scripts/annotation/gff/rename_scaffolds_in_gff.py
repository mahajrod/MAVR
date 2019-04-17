#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import AnnotationsRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_gff", action="store", dest="input_gff", required=True,
                    help="Input GFF file")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-s", "--syn_file", action="store", dest="syn_file", required=True,
                    help="File with scaffold synonyms")

args = parser.parse_args()

AnnotationsRoutines.rename_scaffolds_in_gff(args.input_gff, args.syn_file, args.output_prefix)
