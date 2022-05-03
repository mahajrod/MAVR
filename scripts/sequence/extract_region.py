#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import sys
import argparse

import pandas as pd

from RouToolPa.Parsers.Sequence import CollectionSequence

BED_SCAFFOLD_COL = 0
BED_START_COL = 1
BED_END_COL = 2
BED_ID_COL = 3

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_fasta", action="store", dest="input_fasta", required=True,
                    type=lambda s: s.split(","),
                    help="Input fasta file")
parser.add_argument("-o", "--output", action="store", dest="output", default="stdout",
                    help="Output fasta file with sequences of region")
parser.add_argument("-b", "--bed", action="store", dest="bed",
                    help="BED file with regions")
args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

args.bed = pd.read_csv(args.bed, sep="\t", header=None, index_col=0)

fasta_collection = CollectionSequence(in_file=args.input_fasta, format="fasta")
output_collection = CollectionSequence()
for region_tuple in args.bed.itertuples():
    output_collection.records[str(region_tuple[BED_ID_COL])] = fasta_collection.records[region_tuple[BED_SCAFFOLD_COL]][region_tuple[BED_START_COL]:region_tuple[BED_END_COL]]

output_collection.write(args.output)



