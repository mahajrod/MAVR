#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.Routines import FileRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input GFF file. Default: stdin")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output GFF file. Default: stdout ")

args = parser.parse_args()

with FileRoutines.metaopen(args.input, "r") as in_fd, FileRoutines.metaopen(args.output, "r") as out_fd:
    for line in in_fd:
        if line[0] == "#":
            out_fd.write(line)
            continue
        line_list = line.split("\t")
        if len(line_list) == 1:
            out_fd.write(line)
            continue
        line_list[0] = ".".join(line_list[0].split(".")[:-1])
        out_fd.write(line)

