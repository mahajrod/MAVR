#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.Routines import AnnotationsRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-g", "--gff", action="store", dest="gff", required=True,
                    help="Gff file")
parser.add_argument("-f", "--features", action="store", dest="features",
                    default=[],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of features to count per scaffold. "
                         "If not set all features will be counted")
parser.add_argument("-o", "--output", action="store", dest="output",
                    help="Output file with counts of features. Default: stdout")

args = parser.parse_args()

if args.output is None:
    args.output = sys.stdout

AnnotationsRoutines.count_per_scaffold_feature_number(args.gff, out_file=args.output, feature_type_list=args.features)
