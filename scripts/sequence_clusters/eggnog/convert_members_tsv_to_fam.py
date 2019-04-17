#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.Routines import EggNOGRoutines



parser = argparse.ArgumentParser()

parser.add_argument("-i", "--members_file", action="store", dest="members_file", required=True,
                    help="Members .tsv file from EggNOG")
parser.add_argument("-o", "--output", action="store", dest="output", default="stdout",
                    help="File to write output")

args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

EggNOGRoutines.convert_members_tsv_to_fam(args.members_file, args.output)
