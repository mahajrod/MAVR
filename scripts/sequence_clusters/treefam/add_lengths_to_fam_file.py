#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.Routines import TreeFamRoutines



parser = argparse.ArgumentParser()

parser.add_argument("-i", "--fam_file", action="store", dest="fam_file", required=True,
                    help="File with families")
parser.add_argument("-l", "--length_file", action="store", dest="length_file", required=True,
                    help="File with lengths of members")
parser.add_argument("-o", "--output", action="store", dest="output", default="stdout",
                    help="File to write output")

args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

TreeFamRoutines.add_length_to_fam_file(args.fam_file, args.length_file, out_fd, close_after_if_file_object=True)
