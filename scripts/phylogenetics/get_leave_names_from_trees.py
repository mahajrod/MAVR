#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import argparse
from ete2 import Tree
from RouToolPa.Routines.File import check_path, split_filename
from RouToolPa.Collections.General import SynDict



parser = argparse.ArgumentParser()

parser.add_argument("-i", "--tree_dir", action="store", dest="tree_dir", required=True, type=check_path,
                    help="Directory with trees")
parser.add_argument("-f", "--tree_format", action="store", dest="tree_format", default=1, type=int,
                    help="Format of input trees")
parser.add_argument("-o", "--output_file", action="store", dest="output_file", default="stdout",
                    help="Output file with leaves of trees. Default: stdout")

args = parser.parse_args()

out_fd = sys.stdout if args.output_file == "stdout" else open(args.output_file, "w")

tree_files_list = os.listdir(args.tree_dir)

names_dict = SynDict()

for tree_file in tree_files_list:
    tree_name = split_filename(tree_file)[1]
    with open("%s%s" % (args.tree_dir, tree_file), "r") as tree_fd:
        tree = Tree(tree_fd.readline().strip(), format=args.tree_format)
    leaves_list = []
    for node in tree.traverse():
        if node.is_leaf():
            leaves_list.append(node.name)
    names_dict[tree_name] = leaves_list

names_dict.write(args.outp_fd, splited_values=True)
if args.output_file != "stdout":
    out_fd.close()