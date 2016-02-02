#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

from Routines import AlignmentRoutines, File


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    type=lambda s: File.make_list_of_path_to_files(s.split(",")),
                    help="Comma-separated list of files/directories with alignments")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="File to write merged alignment")
parser.add_argument("-c", "--coordinates_file", action="store", dest="coords_file",
                    help="File to write file with coordinates of alignments in merged alignment")

args = parser.parse_args()

AlignmentRoutines.merge_alignment(args.input, args.output, args.coords_file)
