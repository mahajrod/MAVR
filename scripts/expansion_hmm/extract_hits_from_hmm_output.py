#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import argparse

from Bio import SearchIO

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True, type=lambda s: s.split(","),
                    help="Comma-separated list of files with hmm output")
parser.add_argument("-f", "--format", action="store", dest="format", required=True,
                    help="Format of input hmm file.")
parser.add_argument("-o", "--output_file", action="store", dest="output", default="stdout",
                    help="Output file")
args = parser.parse_args()

hmm_dict = SearchIO.index_db("temp.idx", args.input, args.format)
out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")
out_fd.write("#query\thit\tevalue\tbitscore\n")
for family in hmm_dict:
    #print hmm_dict[key]
    for hit in hmm_dict[family]:
        if hit.is_included:
            out_fd.write("%s\t%s\t%s\t%s\n" % (family, hit.id, hit.evalue, hit.bitscore))
if args.output != "stdout":
    out_fd.close()

os.remove("temp.idx")