#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import AnnotationsRoutines
from RouToolPa.Collections.General import IdList


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_gff", action="store", dest="input_gff", required=True,
                    help="Input .gff file")
parser.add_argument("-f", "--value_file", action="store", dest="value_file", required=True,
                    help="Value with values to seek for")
parser.add_argument("-o", "--output_gff", action="store", dest="output_gff", required=True,
                    help="Output .gff file")
parser.add_argument("-d", "--description_fields", action="store",
                    dest="field_id_list",
                    type=lambda s: s.split(","), required=True,
                    help="Comma-separated list of fields in gff description to check")



args = parser.parse_args()

value_list = IdList(filename=args.value_file)
AnnotationsRoutines.extract_gff_records_by_description_value(args.input_gff, args.output_gff, args.field_id_list, value_list,
                                                retain_comments=False)
