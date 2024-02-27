#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import bz2
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
                    help="Input bed12 file with query_annotations (query_annotations.bed). Default: stdin")
parser.add_argument("-d", "--transcript_ids", action="store", dest="transcript_ids", required=True,
                    help="File with transcript ids extracted from query_isoforms.tsv")
parser.add_argument("-f", "--field_idx", action="store", dest="field_idx", default=3, type=int,
                    help="Index (0-based) of field to use for filtering. Default: 3")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output filtered bed12 file. Default: stdout")

args = parser.parse_args()

transcript_set = set(pd.read_csv(args.transcript_ids, sep="\t", header=None).squeeze("columns"))
#print(pd.read_csv(args.transcript_ids, sep="\t", header=None))
#print(transcript_set)
with metaopen(args.input, "r") as in_fd, metaopen(args.output, "w") as out_fd:
    for line in in_fd:
        if line.split("\t")[args.field_idx] in transcript_set:
            out_fd.write(line)

