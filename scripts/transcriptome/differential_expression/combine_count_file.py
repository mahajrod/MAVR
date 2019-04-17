#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Expression import HTSeq


parser = argparse.ArgumentParser()

parser.add_argument("-f", "--file_list", action="store", dest="file_list", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of files with counts")
parser.add_argument("-o", "--output_file", action="store", dest="output", required=True,
                    help="Output file with combined counts")
parser.add_argument("-s", "--sample_list", action="store", dest="sample_list",
                    type=lambda s: s.split(","),
                    help="Comma-separated list of sample names."
                         "If not set filenames will be treated as sample names")


args = parser.parse_args()

HTSeq.combine_count_files(args.file_list, args.output, sample_name_list=args.sample_list)
