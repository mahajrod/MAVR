#!/usr/bin/env python
__author__ = 'mahajrod'

import bz2
import sys
import gzip
import argparse

import pandas as pd

if sys.version_info[0] == 3:
    from io import TextIOWrapper as file


def metaopen(filename, flags, buffering=None, compresslevel=5):
    if not isinstance(filename, str): # or isinstance(filename, gzip.GzipFile) or isinstance(filename, bz2.BZ2File):
        if isinstance(filename, file):
            return filename
        else:
            raise ValueError("ERROR!!! Not file object or str: {}".format(str(filename)))
    elif filename[-3:] == ".gz":
        return gzip.open(filename, flags + ("t" if "b" not in flags else ""), compresslevel=compresslevel)
    elif filename[-4:] == ".bz2":
        return bz2.open(filename, flags + ("t" if "b" not in flags else ""), compresslevel=compresslevel)
    else:
        if buffering is not None:
            return open(filename, flags, buffering=buffering)
        else:
            return open(filename, flags)

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input sam file. Default: stdin")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output file. Default: stdout")
parser.add_argument("-d", "--read_id_file", action="store", dest="output", default=None,
                    help="File with read ids to extract. Default: not set, i.e all the reads will be included ")

args = parser.parse_args()

with metaopen(args.input, "r") as in_fd, metaopen(args.output, "w") as out_fd:
    if args.read_id_file is None:
        for line in in_fd:
            out_fd.write(line)
    else:
        read_series = pd.read_csv(args.read_id_file, sep="\t", header=None).squeeze("columns")
        for line in in_fd:
            if line[0] == "@":
                out_fd.write(line)
                continue
            if line.split("\t")[0] in read_series:
                out_fd.write(line)
