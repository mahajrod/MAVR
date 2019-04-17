#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Parsers.PAML import CodeMLReport

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
codeml_report.extract_trees(args.out_prefix)
