#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.GeneralRoutines import FileRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input fasta file. Default: stdin")
parser.add_argument("-s", "--separator", action="store", dest="separator", default=".",
                    help="Separator between id and version in the scaffold id. "
                         "Everything in the scaffold id after last separator will be trimmed. "
                         "Description will be retained. "
                         "However, TABs and other spacing symbols will be replaced by spaces. Default: '.'")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output fasta file with removed versions of scaffolds. Default: stdout")

args = parser.parse_args()

with FileRoutines.metaopen(args.input, "r") as in_fd, \
     FileRoutines.metaopen(args.output, "w") as out_fd:
    for line in in_fd:
        if line[0] != ">":
            out_fd.write(line)
        else:
            line_list = line.strip().split()
            line_list[0] = line_list[0].rsplit(args.separator, 1)[0]
            out_fd.write(" ".join(line_list) + "\n")

