#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from ete2 import Tree
from RouToolPa.Routines.File import read_synonyms_dict

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_tree_file", action="store", dest="input_tree_file",
                    help="File with input trees")
parser.add_argument("-o", "--output_tree_file", action="store", dest="output_tree_file",
                    help="File with output trees")
parser.add_argument("-f", "--tree_format", action="store", dest="tree_format", type=int, default=0,
                    help="""Format of input trees. Allowed formats:
0 	flexible with support values (default)
1 	flexible with internal node names
2 	all branches + leaf names + internal supports
3 	all branches + all names
4 	leaf branches + leaf names
5 	internal and leaf branches + leaf names
6 	internal branches + leaf names
7 	leaf branches + all names
8 	all names
9 	leaf names
100 	topology only""")

parser.add_argument("-s", "--file_with_synonyms", action="store", dest="file_with_synonyms", required=True,
                    help="File with synonyms. ")

args = parser.parse_args()

out_fd = open(args.output_tree_file, "w")
tree_index = 1

synonyms_dict = read_synonyms_dict(args.file_with_synonyms, separator="\t", header=False)

with open(args.input_tree_file, "r") as in_fd:
    for line in in_fd:
        children_index = 1
        tree_line = line.strip()
        tree = Tree(tree_line, format=args.tree_format)

        print("Handling tree %i" % tree_index)

        for node in tree.traverse():
            if node.name in synonyms_dict:
                node.name = synonyms_dict[node.name]

        out_fd.write(tree.write(format=args.tree_format)) #, features=tree.features - set(["support", "name"])))

        tree_index += 1

out_fd.close()