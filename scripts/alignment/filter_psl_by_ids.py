#!/usr/bin/env python
__author__ = 'mahajrod'
import sys
import argparse
from RouToolPa.Routines import AlignmentRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input PSL file with. Default: stdin")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output file. Default: stdout")

parser.add_argument("-w", "--white_query_id_file", action="store", dest="white_query_id_file",
                    help="File with query ids from white list")
parser.add_argument("-b", "--black_query_id_file", action="store", dest="black_query_id_file",
                    help="File with query ids from black list")

parser.add_argument("-x", "--white_target_id_file", action="store", dest="white_target_id_file",
                    help="File with target ids from white list")
parser.add_argument("-y", "--black_target_id_file", action="store", dest="black_target_id_file",
                    help="File with target ids from black list")


#parser.add_argument("-n", "--no_header", action="store_true", dest="no_header",
#                    help="Don't add header to output")

args = parser.parse_args()

AlignmentRoutines.filter_psl_by_ids_from_file(args.input, args.output,
                                              white_query_id_file=args.white_query_id_file,
                                              black_query_id_file=args.black_query_id_file,
                                              white_target_id_file=args.white_target_id_file,
                                              black_target_id_file=args.black_target_id_file)
