#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import AnnotationsRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_gff", action="store", dest="input_gff", required=True,
                    help="Input .gff file")
parser.add_argument("-o", "--output_gff", action="store", dest="output_gff", required=True,
                    help="Output .gff file")
parser.add_argument("-s", "--syn_file_file", action="store", dest="syn_file", required=True,
                    help="File with synonyms of region names")

args = parser.parse_args()

AnnotationsRoutines.replace_region_names_in_gff(args.input_gff, args.syn_file, args.output_gff)
