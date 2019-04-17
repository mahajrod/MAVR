#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import AnnotationsRoutines

parser = argparse.ArgumentParser()


parser.add_argument("-i", "--input_gff", action="store", dest="input_gff", required=True,
                    help="input gff file")
parser.add_argument("-o", "--output_gff", action="store", dest="output_gff", required=True,
                    help="Output fixed gff file")
parser.add_argument("-f", "--feature_type", action="store", dest="feature_type", required=True,
                    help="Feature type to use")
args = parser.parse_args()

AnnotationsRoutines.fix_absent_feature_type_field(args.input_gff, args.output_gff, args.feature_type)
