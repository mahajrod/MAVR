#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import argparse

from Bio import SearchIO

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input file with hmm output")
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
    try:
        if hmm_dict[family][0].is_included:
            out_fd.write("%s\t%s\t%s\t%s\n" % (family, hmm_dict[family][0].id, hmm_dict[family][0].evalue,
                                               hmm_dict[family][0].bitscore))
    except:
        print hmm_dict[family]    

if args.output != "stdout":
    out_fd.close()

os.remove("temp.idx")