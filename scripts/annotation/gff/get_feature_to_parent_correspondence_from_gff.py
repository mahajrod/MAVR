#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import AnnotationsRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_gff", action="store", dest="input_gff", required=True,
                    help="Input .gff file")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="File to write output")
parser.add_argument("-t", "--feature_types", action="store", dest="feature_types",
                    type=lambda s: s.split(","), default=["mRNA", ],
                    help="Comma-separated list of feature types to count. "
                         "Default: mRNA")
parser.add_argument("-d", "--id_entry", action="store", dest="id_entry", default="ID",
                    help="Id entry. Default: ID")
parser.add_argument("-p", "--parental_id_entry", action="store", dest="parental_id_entry", default="Parent",
                    help="Parental id entry. Default: Parent")

args = parser.parse_args()

AnnotationsRoutines.get_feature_to_parent_correspondence_from_gff(args.input_gff, args.output,
                                                                  feature_list=args.feature_types,
                                                                  id_entry=args.id_entry,
                                                                  parental_id_entry=args.parental_id_entry)
