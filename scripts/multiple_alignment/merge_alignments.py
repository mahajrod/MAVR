#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import MultipleAlignmentRoutines, FileRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    type=lambda s: FileRoutines.make_list_of_path_to_files(s.split(",")),
                    help="Comma-separated list of files/directories with alignments")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="File to write merged alignment")
parser.add_argument("-c", "--coordinates_file", action="store", dest="coords_file", required=True,
                    help="File to write file with coordinates of alignments in merged alignment")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of alignments")

args = parser.parse_args()

MultipleAlignmentRoutines.merge_alignment(args.input, args.output, args.coords_file, format=args.format)
