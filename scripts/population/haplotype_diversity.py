#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument("-c", "--counts", action="store", dest="counts", required=True,
                    help="Tab-separated file with counts")
parser.add_argument("-p", "--output_prefix", action="store", dest="output_prefix",
                    help="Prefix for output_files")

args = parser.parse_args()

counts_df = pd.read_csv(args.counts, names=("haplotype", 'counts'), converters={"counts": int}, sep="\t")
total_samples = counts_df["counts"].sum()
counts_df["frequency"] = counts_df["counts"] / total_samples
counts_df["frequence2"] = counts_df["frequency"] * counts_df["frequency"]

H = float(total_samples) / float(total_samples - 1) * (1 - counts_df["frequence2"].sum())

print(H)

if args.output_prefix:
    counts_df.to_csv("%s.df" % args.output_prefix, sep="\t")
