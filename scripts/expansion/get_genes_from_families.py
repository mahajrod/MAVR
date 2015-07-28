#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse

from CustomCollections.GeneralCollections import IdList
from Routines.File import read_synonyms_dict

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input file with families")
parser.add_argument("-id", "--id_file", action="store", dest="id_file", required=True,
                    help="File with ids of families.")
parser.add_argument("-o", "--output_file", action="store", dest="output",
                    help="Output file")
args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

id_list = IdList
id_list = IdList.read(args.id_file)

families = read_synonyms_dict(args.input,separator="\t", split_values=True, values_separator=",")

with open(args.output, "w") as out_fd:
    for fam_id in id_list:
        for gene_id in families[fam_id]:
            out_fd.write(gene_id + "\n")
if args.output != "output":
    out_fd.close()
