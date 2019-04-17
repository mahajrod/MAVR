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
parser.add_argument("-d", "--white_id_file", action="store", dest="white_id_file",
                    help="File with ids from white list. If set other ids are ignored")
parser.add_argument("-m", "--max_hits_per_query", action="store", dest="max_hits_per_query", type=int,
                    help="Maximum hits per query")
args = parser.parse_args()

top_hits_gff = "%s.target.top_hits.gff" % args.output_prefix
secondary_hits_gff = "%s.target.secondary_hits.gff" % args.output_prefix

Exonerate.extract_top_hits_from_target_gff(args.input, top_hits_gff, secondary_hits_gff,
                                           id_white_list_file=args.white_id_file,
                                           max_hits_per_query=args.max_hits_per_query)
