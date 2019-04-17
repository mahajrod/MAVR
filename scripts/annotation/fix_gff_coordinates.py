#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import AnnotationsRoutines

parser = argparse.ArgumentParser()


parser.add_argument("-i", "--input_gff", action="store", dest="input_gff", required=True,
                    help="input gff file")
parser.add_argument("-o", "--output_gff", action="store", dest="output_gff", required=True,
                    help="Output fixed gff file")

args = parser.parse_args()

AnnotationsRoutines.fix_gff_coordinates_order(args.input_gff, args.output_gff)
