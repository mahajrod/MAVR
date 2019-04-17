#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import AnnotationsRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_gff", action="store", dest="input_gff", required=True,
                    help="Input .gff file")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Output .gff file")

parser.add_argument("-t", "--feature_types", action="store", dest="feature_types",
                    type=lambda s: s.split(","), default=["CDS"],
                    help="Comma-separated list of feature types to extract. "
                         "Default: CDS only")

parser.add_argument("-u", "--unification_key", action="store",
                    dest="unification_key",
                    default="Parent",
                    help="Annotation entry to use for unification. Default: Parent")

args = parser.parse_args()

AnnotationsRoutines.get_feature_dict(args.input_gff,
                                     output_prefix=args.output_prefix,
                                     feature_type_list=args.feature_types,
                                     unification_key=args.unification_key)
