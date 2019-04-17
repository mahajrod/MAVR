#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import MultipleAlignmentRoutines, FileRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    type=lambda s: FileRoutines.make_list_of_path_to_files(s.split(",")),
                    help="Comma-separated list of files/directories with alignments")
parser.add_argument("-o", "--output_directory", action="store", dest="output_dir", default="./",
                    help="Output directory to write resulting files. Default - current directory")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of alignments")
parser.add_argument("-g", "--gap_symbol", action="store", dest="gap_symbol", default="-",
                    help="Gap symbol. Default - '-'")

args = parser.parse_args()

for alignment_file in args.input:
    alignment_name_list = FileRoutines.split_filename(alignment_file)
    output_file = "%s/%s.position_matrix" % (args.output_dir, alignment_name_list[1])
    MultipleAlignmentRoutines.get_position_presence_matrix_fom_file(alignment_file, output_file,
                                                                    format=args.format, gap_symbol=args.gap_symbol)
