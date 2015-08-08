#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import sys
import argparse

from ete2 import Tree

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_tree_file", action="store", dest="input_tree_file",
                    help="File with input tree")
parser.add_argument("-f", "--tree_format", action="store", dest="tree_format", default=1, type=int,
                    help="Format of input tree")
parser.add_argument("-o", "--output_file", action="store", dest="output_file", default="stdout",
                    help="Output file with names of nodes. Default: stdout")
parser.add_argument("-m", "--extraction_mode", action="store", dest="extraction_mode", default="leaves",
                    help="Mode used to extract name: all, leaves(default), internal")

args = parser.parse_args()

out_fd = sys.stdout if args.output_file == "stdout" else open(args.output_file, "w")
tree_index = 1
if args.extraction_mode == "all":
    extraction_func = lambda node: True
elif args.extraction_mode == "leaves":
    extraction_func = lambda node: node.is_leaf()
elif args.extraction_mode == "internal":
    extraction_func = lambda node: not node.is_leaf()


with open(args.input_tree_file, "r") as in_fd:
    for line in in_fd:
        tree_line = line.strip()
        tree = Tree(tree_line, format=args.tree_format)
        out_fd.write("#Tree %i\n" % tree_index)
        for node in tree.traverse():
            if extraction_func(node):
                out_fd.write("%s\n" % node.name)
        tree_index += 1

if args.output_file != "stdout":
    out_fd.close()