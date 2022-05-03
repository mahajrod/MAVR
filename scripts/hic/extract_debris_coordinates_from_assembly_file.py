#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.Collections.General import SynDict


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input .assembly file")

parser.add_argument("-o", "--output", action="store", dest="output", default="stdout",
                    help="Output .bed file with coordinates od debris fragments")

args = parser.parse_args()

with open(args.input, "r") as in_fd, open(args.output, "w") as out_fd:
    prev_id = ""
    prev_end = 0
    for line in in_fd:
        if line[0] != ">": # ignore assembly lines
            continue
        line_list = line.strip()[1:].split(" ")
        current_id_list = line_list[0].split(":::")

        length = int(line_list[-1])
        if current_id_list[0] != prev_id:
            prev_id = current_id_list[0]
            prev_end = 0

        if current_id_list[-1] == "debris":
            out_fd.write("{0}\t{1}\t{2}\t{3}\n".format(current_id_list[0], prev_end, prev_end + length,
                                                       line_list[-2]))

        prev_end += length
