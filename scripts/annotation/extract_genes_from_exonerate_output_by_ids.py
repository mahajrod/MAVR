#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Annotation import Exonerate
from RouToolPa.Routines.File import make_list_of_path_to_files


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    type=lambda s: make_list_of_path_to_files(s.split(",")),
                    help="Input comma-separated list of files/directories with exonerate output")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    help="Prefix of output files")
parser.add_argument("-d", "--id_file", action="store", dest="id_file",
                    help="File with ids of sequences to extract")
args = parser.parse_args()

extracted_gff = "%s.filtered.gff" % args.output_prefix
filtered_out_gff = "%s.filtered_out.gff" % args.output_prefix

Exonerate.extract_annotation_by_refence_id(args.input, args.id_file, extracted_gff, filtered_out_gff)
