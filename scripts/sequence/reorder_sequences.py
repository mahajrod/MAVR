#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse

import pandas as pd

from RouToolPa.Parsers.Sequence import CollectionSequence


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input fasta file. May be gziped. Default: stdin")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output reordered fasta file. Default: stdout")
parser.add_argument("-r", "--orderlist", action="store", dest="orderlist", default=None,
                    type=lambda s: pd.read_csv(s, sep="\t", header=None).squeeze(),
                    help="File with order list of the sequences (one id per line). Default: not set")
parser.add_argument("-b", "--by", action="store", dest="by", default="length",
                    help="How to sort sequenced. Allowed: 'length' (default), orderlist")



args = parser.parse_args()

sequence_col = CollectionSequence(in_file=args.input)
sequence_col.reorder_records(by=args.by, orderlist=args.orderlist, in_place=True)

sequence_col.write(args.output)

