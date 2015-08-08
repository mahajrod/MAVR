#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse

from CustomCollections.GeneralCollections import IdList
from Routines.File import read_synonyms_dict

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input file with families")
parser.add_argument("-d", "--id_file", action="store", dest="id_file", default=None,
                    help="File with ids of families. If absent genes from all families will be extracted(default).")
parser.add_argument("-o", "--output_file", action="store", dest="output",
                    help="Output file")
args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

families = read_synonyms_dict(args.input, separator="\t", split_values=True, values_separator=",")
if args.id_file:
    id_list = IdList()
    id_list = id_list.read(args.id_file)

with open(args.output, "w") as out_fd:
    for fam_id in id_list if args.id_file else families:
        for gene_id in families[fam_id]:
            out_fd.write(gene_id + "\n")
if args.output != "output":
    out_fd.close()
