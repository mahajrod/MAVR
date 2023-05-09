#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse

from RouToolPa.GeneralRoutines import FileRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input fastq file. Default: stdin")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output file. Default: stdout")

args = parser.parse_args()

with FileRoutines.metaopen(args.input, "r") as in_fd, FileRoutines.metaopen(args.output, "w") as out_fd:
    for line in in_fd:
        read_id = line.strip().split()[0][1:]
        length = len(in_fd.readline().strip())
        out_fd.write("{0}\t{1}\n".format(read_id, length))
        in_fd.readline()
        in_fd.readline()
