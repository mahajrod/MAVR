#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import AnnotationsRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_gff", action="store", dest="input_gff",
                    help="Input .gff file")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    help="Output prefix")
parser.add_argument("-t", "--feature_types", action="store", dest="feature_types",
                    type=lambda s: s.split(","), default=None,
                    help="Comma-separated list of feature types to add aliases. "
                         "Default: all feature types")

args = parser.parse_args()

AnnotationsRoutines.count_total_feature_length_from_gff(args.input_gff, args.output_prefix,
                                                        features_to_count=args.feature_types)
