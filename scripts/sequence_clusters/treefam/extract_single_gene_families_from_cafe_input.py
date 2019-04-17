#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input cafe file")

parser.add_argument("-o", "--output", action="store", dest="output", default="stdout",
                    help="File to write single gene families. Default - stdout")

args = parser.parse_args()

out_file = sys.stdout if args.output == "stdout" else open(args.output, "w")

with open(args.input, "r") as in_fd:
    header_list = in_fd.readline().strip().split("\t")
    out_file.write("\t".join(header_list[1:]) + "\n")
    for line in in_fd:
        line_list = line.strip().split("\t")
        fam_name = line_list[1]
        counts_list = map(int, line_list[2:])
        for count in counts_list:
            if count != 1:
                break
        else:
            out_file.write("%s\t%s\n" % (fam_name, "\t".join(map(str, counts_list))))


