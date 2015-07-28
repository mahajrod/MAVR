#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input file with hclust output")
parser.add_argument("-o", "--output_file", action="store", dest="output",
                    help="Output file")
args = parser.parse_args()

with open(args.input, "r") as in_fd:
    with open(args.output, "w") as out_fd:
        for line in in_fd:
            line_list = line.strip().split("\t")
            out_fd.write("%s\t%s\n" % (line_list[0], line_list[-1][:-1]))


