#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse

from Parsers.PAML import CodeMLReport

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with CodeML report")
parser.add_argument("-o", "--output_file_prefix", action="store", dest="out_prefix", required=True,
                    default="codeml_trees",
                    help="Prefix of output files with trees")
parser.add_argument("-t", "--tree", action="store", dest="tree",
                    help="File with phylogenetic tree")

args = parser.parse_args()

codeml_report = CodeMLReport(args.input, treefile=args.tree)
codeml_report.write_trees(args.out_prefix)
codeml_report.get_feature_values(mode="all")
codeml_report.get_feature_values(mode="leaves")
codeml_report.get_feature_values(mode="internal")
codeml_report.get_all_values("dN_dS_W.t")
codeml_report.get_leaf_values()
if codeml_report.branches_with_positive_selection():
    sys.stderr.write("Presence of branches with positive selection\n")

codeml_report.convert_trees_to_tsv(args.out_prefix)




