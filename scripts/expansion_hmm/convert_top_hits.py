#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import argparse


from collections import OrderedDict

from CustomCollections.GeneralCollections import SynDict
from Routines.File import read_synonyms_dict

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input file with top hits")
parser.add_argument("-e", "--header", action="store_true", dest="header",
                    help="Header is present in input file")
parser.add_argument("-o", "--output_file", action="store", dest="output", default="stdout",
                    help="Output file")
args = parser.parse_args()

syn_dict = read_synonyms_dict(args.input, header=args.header, split_values=False)

hit_dict = SynDict()
for query in syn_dict:
    hit = syn_dict[query]
    if hit not in hit_dict:
        hit_dict[hit] = [query]
    else:
        hit_dict[hit].append(query)

SynDict.write(args.output,splited_values=True)