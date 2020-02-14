#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import AnnotationsRoutines
from RouToolPa.Collections.General import IdList


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_gff", action="store", dest="input_gff", required=True,
                    help="Input .gff file")
parser.add_argument("-a", "--annotation_types", action="store", dest="annotation_types", default=["gene"],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of annotation types to extract "
                         "(including sub annotations). Default: 'gene'")
parser.add_argument("-g", "--feature_ids", action="store", dest="feature_ids",
                    help="Comma-separated list of IDs for features to be extracted")
parser.add_argument("-f", "--feature_id_file", action="store", dest="feature_id_file",
                    help="File with ids of feature to be extracted. Ignored if -g/--feature_ids option is set")
parser.add_argument("-o", "--output_gff", action="store", dest="output_gff", required=True,
                    help="Output .gff file")


args = parser.parse_args()

feature_ids = args.feature_ids.split(",") if args.feature_ids else IdList(filename=args.feature_id_file)

AnnotationsRoutines.extract_annotation_from_gff(args.input_gff,
                                                args.feature_ids,
                                                args.annotation_types,
                                                args.output_gff)
