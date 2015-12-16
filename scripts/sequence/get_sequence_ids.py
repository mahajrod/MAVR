#!/usr/bin/env python
__author__ = 'mahajrod'
import os
import sys
import argparse

from Bio import SeqIO


parser = argparse.ArgumentParser()

parser.add_argument("-f", "--format", action="store", dest="format",
                    default="fasta", help="format of file with sequences - default: fasta.")
parser.add_argument("-i", "--input", action="store", dest="input",
                    help="File with sequences")
parser.add_argument("-o", "--output", action="store", dest="output",
                    default="stdout", help="output file - default: stdout")

args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

for record in SeqIO.parse(args.input, format=args.format):
    out_fd.write("%s\n" % record.id)

if args.output != "stdout":
    out_fd.close()
