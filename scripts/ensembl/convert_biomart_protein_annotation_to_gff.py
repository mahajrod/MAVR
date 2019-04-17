#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import EnsemblRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with protein annotations generated from Biomart")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output gff file")
parser.add_argument("-s", "--separator", action="store", dest="separator", default="\t",
                    help="Column separator in input file. Default: '\t'")
parser.add_argument("-e", "--extraction_mode", action="store", dest="extraction_mode", default="pfam",
                    help="Type of annotation to extract. Allowed: pfam(default).")

args = parser.parse_args()

EnsemblRoutines.convert_biomart_protein_annotation_to_gff(args.input, args.output, separator=args.separator,
                                                          extraction_mode=args.extraction_mode)
