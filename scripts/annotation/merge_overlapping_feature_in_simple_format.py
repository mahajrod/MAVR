#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import AnnotationsRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    type=AnnotationsRoutines.make_list_of_path_to_files_from_string,
                    help="Comma-separated list of files or directories with files"
                         "containing annotations with three columnst: <scaffold>, <start>, <stop>")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file")

parser.add_argument("-r", "--scaffold_column_id", action="store", dest="scaffold_column_id",
                    type=int, default=0,
                    help="0-based index of column with scaffold id. Default: 0")
parser.add_argument("-s", "--start_column_id", action="store", dest="start_column_id",
                    type=int, default=1,
                    help="0-based index of column with feature start. Default: 1")
parser.add_argument("-e", "--end_column_id", action="store", dest="end_column_id",
                    type=int, default=2,
                    help="0-based index of column with feature end. Default: 2")
parser.add_argument("-n", "--coordinates_type", action="store", dest="coordinates_type",
                    default="1-based",
                    help="Type of coordinates. Allowed: 0-based, 1-based(default)")

args = parser.parse_args()

AnnotationsRoutines.merge_overlapping_feature_in_simple_format(args.input,
                                                               args.scaffold_column_id,
                                                               args.start_column_id,
                                                               args.end_column_id,
                                                               output_file=args.output,
                                                               output_separator="\t",
                                                               comments_prefix="#",
                                                               input_separator="\t",
                                                               coordinates_type=args.coordinates_type,
                                                               return_seqfeature_dict=False,
                                                               feature_type=None)
