#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.Routines import EggNOGRoutines



parser = argparse.ArgumentParser()

parser.add_argument("-i", "--fam_file", action="store", dest="fam_file", required=True,
                    help="File with families obtained from hits to hmm profiles of EggNOG")
parser.add_argument("-o", "--output", action="store", dest="output", default="stdout",
                    help="File to write output")

args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

EggNOGRoutines.edit_profile_names_in_fam_file(args.fam_file, args.output)
