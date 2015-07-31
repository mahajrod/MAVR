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
parser.add_argument("-n", "--not_found", action="store", dest="not_found", default="not_found.ids",
                    help="File to write ids of not found queries")
parser.add_argument("-g", "--not_significant", action="store", dest="not_significant", default="not_significant.hits",
                    help="File to write ids of quaries with not significant hits")
args = parser.parse_args()

hmm_dict = SearchIO.index_db("temp.idx", args.input, args.format)
out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")
out_fd.write("#query\thit\tevalue\tbitscore\n")

nf_fd = open(args.not_found, "w")
ns_fd = open(args.not_significant, "w")
for query in hmm_dict:
    #print hmm_dict[key]
    print hmm_dict[query].hits
    if hmm_dict[query]:
        if hmm_dict[query][0].is_included:
            out_fd.write("%s\t%s\t%s\t%s\n" % (query, hmm_dict[query][0].id, hmm_dict[query][0].evalue,
                                               hmm_dict[query][0].bitscore))
        else:
            ns_fd.write("%s\t%s\t%s\t%s\n" % (query, hmm_dict[query][0].id, hmm_dict[query][0].evalue,
                                               hmm_dict[query][0].bitscore))
    else:
        out_fd.write("%s\n" % query)

if args.output != "stdout":
    out_fd.close()

os.remove("temp.idx")

nf_fd.close()
ns_fd.close