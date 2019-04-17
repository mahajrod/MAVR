#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from random import choice
from RouToolPa.Collections.General import IdList, SynDict

# from RouToolPa.Routines.File import read_synonyms_dict

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input file with node information")
parser.add_argument("-r", "--reference_fam", action="store", dest="reference_fam", required=True,
                    help="Reference family file")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", default="stdout",
                    help="Prefix of output file")
args = parser.parse_args()

out_fd = sys.stdout if args.output_prefix == "stdout" else open("%s_reference_random_genes.ids" % args.output_prefix, "w")

reference_families = SynDict()
reference_families.read(args.reference_fam, separator="\t", split_values=True, values_separator=",")

node_family_ids = IdList()
node_family_ids.read(args.input, header=True, column_number=0, column_separator="\t")

reference_random_genes = SynDict()

for family_id in node_family_ids:
    if family_id not in reference_families:
        reference_random_genes[family_id] = "."
    else:
        reference_random_genes[family_id] = choice(reference_families[family_id])

reference_random_genes.write("%s_reference_random_genes.t" % args.output_prefix)

for family_id in reference_random_genes:
    if reference_random_genes[family_id] != ".":
        out_fd.write("%s\n" % reference_random_genes[family_id])