#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse

from Routines import AnnotationsRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_gff", action="store", dest="input_gff", required=True,
                    help="Input GFF file")
parser.add_argument("-o", "--output", action="store", dest="output_bed", required=True,
                    help="File to write output BED file")
parser.add_argument("-t", "--feature_types", action="store", dest="feature_types",
                    type=lambda s: s.split(","), default=[],
                    help="Comma-separated list of feature types to write in output file "
                         "Default: all")
"""
parser.add_argument("-d", "--id_entry", action="store", dest="id_entry", default="ID",
                    help="Id entry. Default: ID")
parser.add_argument("-p", "--parental_id_entry", action="store", dest="parental_id_entry", default="Parent",
                    help="Parental id entry. Default: Parent")
"""

args = parser.parse_args()

AnnotationsRoutines.convert_gff_to_simple_bed(args.input_gff, args.output_bed, feature_type_list=args.feature_types)
