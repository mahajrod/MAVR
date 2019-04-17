#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from ete2 import Tree
from RouToolPa.Routines.File import read_ids

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_tree_file", action="store", dest="input_tree_file",
                    help="File with input trees to collapse")
parser.add_argument("-o", "--output_tree_file", action="store", dest="output_tree_file",
                    help="File with otput collapsed trees")
parser.add_argument("-f", "--input_tree_format", action="store", dest="input_tree_format", type=int, default=0,
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

parser.add_argument("-d", "--id_file", action="store", dest="id_file",
                    help="File with ids of nodes to remove")

args = parser.parse_args()

id_list = read_ids(args.id_file)

print("Nodes with ids present in %s file will be removed" % args.id_file)

tree_index = 1
with open(args.input_tree_file, "r") as in_fd:
    with open(args.output_tree_file, "w") as out_fd:
        for line in in_fd:
            tree_line = line.strip()
            tree = Tree(tree_line, format=args.input_tree_format)
            print("Totaly %i leaves in tree %i" % (len(tree), tree_index))
            #print(tree.write())
            for node in tree.traverse():
                #if node.is_leaf():
                #    print node.features
                # node.name
                if node.name in id_list:
                    node.delete()
            #print(tree.write())
            print("Totaly %i leaves were retained in tree %i" % (len(tree), tree_index))
            out_fd.write(tree.write(format=args.input_tree_format))
            out_fd.write("\n")
            tree_index += 1
