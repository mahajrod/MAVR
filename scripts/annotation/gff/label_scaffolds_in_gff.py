#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import AnnotationsRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_gff", action="store", dest="input_gff", required=True,
                    help="Input GFF file")
parser.add_argument("-o", "--output_gff", action="store", dest="output_gff", required=True,
                    help="Output gff file")
parser.add_argument("-l", "--label", action="store", dest="label", required=True,
                    help="Label to use")
parser.add_argument("-s", "--separator", action="store", dest="separator", default=".",
                    help="Separator to use. Default: '.'")
args = parser.parse_args()

AnnotationsRoutines.label_scaffolds_in_gff(args.input_gff, args.label, args.output_gff, separator=args.separator)
