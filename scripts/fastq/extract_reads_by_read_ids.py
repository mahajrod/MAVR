#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse

import pandas as pd

from RouToolPa.GeneralRoutines import FileRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdout,
                    help="Input fastq file")
parser.add_argument("-o", "--output", action="store", dest="output", default="stdout",
                    help="Output fastq file")
parser.add_argument("-d", "--id_file", action="store", dest="id_file",
                    help="File with ids of reads to extract")
parser.add_argument("-e", "--excluded", action="store", dest="excluded", default=None,
                    help="File to write excluded reads. If not set corresponding reads will be dropped."
                         " Default: not set ")
#parser.add_argument("-v", "--invert_match", action="store", dest="invert_match", default=False,
#                    help="Invert match, i.e. exclude reads instead of storing. Default: False")
args = parser.parse_args()

read_id_set = set(pd.read_csv(args.input, sep="\t", header=None))
print(len(read_id_set))

if args.excluded is None:
    with FileRoutines.metaopen(args.input, "r") as in_fd, FileRoutines.metaopen(args.output, "w") as out_fd:
        for line in in_fd:
            if in_fd.split()[1:] in read_id_set:
                out_fd.write(line)
                out_fd.write(in_fd.readline())
                out_fd.write(in_fd.readline())
                out_fd.write(in_fd.readline())
            else:
                in_fd.readline()
                in_fd.readline()
                in_fd.readline()
else:
    with FileRoutines.metaopen(args.input, "r") as in_fd, FileRoutines.metaopen(args.output, "w") as out_fd, \
         FileRoutines.metaopen(args.excluded, "w") as excl_fd:
        for line in in_fd:
            if in_fd.split()[1:] in read_id_set:
                out_fd.write(line)
                out_fd.write(in_fd.readline())
                out_fd.write(in_fd.readline())
                out_fd.write(in_fd.readline())
            else:
                excl_fd.write(line)
                excl_fd.write(in_fd.readline())
                excl_fd.write(in_fd.readline())
                excl_fd.write(in_fd.readline())


