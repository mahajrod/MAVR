#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import AnnotationsRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_gff", action="store", dest="input_gff",
                    help="Input .gff file")
parser.add_argument("-s", "--syn_file", action="store", dest="syn_file",
                    help="File with feature synonyms")
parser.add_argument("-o", "--output_gff", action="store", dest="output_gff",
                    help="Output .gff file")
parser.add_argument("-t", "--feature_types", action="store", dest="feature_types",
                    type=lambda s: s.split(","), default=None,
                    help="Comma-separated list of feature types to add aliases. "
                         "Default: all feature types")
parser.add_argument("-n", "--feature_name_fields", action="store",
                    dest="feature_name_fields",
                    type=lambda s: s.split(","), default=["ID", "Name"],
                    help="Comma-separated list of feature name fields in gff description"
                         "Default: ID, Name")
parser.add_argument("-a", "--alias_field", action="store",
                    dest="alias_field",
                    default="Alias",
                    help="Name of field in gff description to add aliases. "
                         "If this field is absent it will be created."
                         "Default: Alias")
parser.add_argument("-k", "--key_column", action="store",
                    dest="key_column", type=int,
                    default=0,
                    help="Key column in synonym file(0-based). Default: 0")
parser.add_argument("-v", "--value_column", action="store",
                    dest="value_column", type=int,
                    default=1,
                    help="Value column in synonym file(0-based). Default: 1")


args = parser.parse_args()

AnnotationsRoutines.add_alias_to_feature(args.input_gff,
                                         args.output_gff,
                                         args.syn_file,
                                         feature_type_list=args.feature_types,
                                         name_field_list=args.feature_name_fields,
                                         alias_field=args.alias_field,
                                         key_column=args.key_column,
                                         value_column=args.value_column)
