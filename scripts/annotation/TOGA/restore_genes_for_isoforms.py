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

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with isoforms ids extracted from query_annotations.bed")
parser.add_argument("-r", "--reference", action="store", dest="reference", required=True,
                    help="Input file with reference isoforms")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file with restored gene to isoform relations.")

args = parser.parse_args()

isoforms = pd.read_csv(args.input, sep="\t", header=None, names=["transcript"])
print(isoforms)
isoforms["reference_transcript"] = isoforms["transcript"].apply(lambda entry: entry.rsplit(".", maxsplit=1)[0])


reference_isoform_df = pd.read_csv(args.reference,
                                   sep="\t",
                                   header=None, names=["gene", "reference_transcript"],
                                   index_col="reference_transcript")

isoform_to_gene_dict = reference_isoform_df.to_dict()["gene"]

isoforms["gene"] = isoforms["reference_transcript"].apply(lambda s: isoform_to_gene_dict[s] if s in isoform_to_gene_dict else pd.NA)

isoforms[["gene", "transcript"]].to_csv(args.output, sep="\t", header=False, index=False)
