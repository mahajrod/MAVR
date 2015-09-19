#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os
import sys
import argparse

from Bio import SeqIO

from Routines.File import make_list_of_path_to_files

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    type=lambda s: make_list_of_path_to_files(s.split(",")),
                    help="Comma-separated list of genbank files/directories with microsatellite data")
parser.add_argument("-o", "--output", action="store", dest="output", default="stdout",
                    help="Output file. Default: stdout")

args = parser.parse_args()

sequence_dict = SeqIO.index_db("temp.idx", args.input, format="genbank")
out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

for record in


if args.output != "output":
    out_fd.close()
