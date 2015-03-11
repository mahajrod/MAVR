#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

from ete2 import Tree


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

parser.add_argument("-s", "--support_threshold", action="store", dest="minimum_support", default=50,
                    help="Minimun support value for nodes to retain")
parser.add_argument("-y", "--support_type", action="store", dest="support_type", default="int",
                    help="Type of support values - 'int' or 'float'. Default - 'int'")

args = parser.parse_args()

threshold_value = int(args.minimum_support) if args.support_type == "int" else float(args.minimum_support) \
    if args.support_type == "float" else None

print("Nodes with support less than %s will be collapsed" % str(threshold_value))
if threshold_value is None:
    raise ValueError("Wrong support type is set")
tree_index = 1
with open(args.input_tree_file, "r") as in_fd:
    with open(args.output_tree_file, "w") as out_fd:
        for line in in_fd:
            tree_line = line.strip()
            tree = Tree(tree_line, format=args.input_tree_format, support=100 if args.support_type == "int" else 1.0)
            print("Totaly %i leaves in tree %i" % (len(tree), tree_index))
            #print(tree.write())
            for node in tree.traverse():
                if node.is_root() or node.is_leaf():
                    #print(node.support)
                    #print(tree.write())
                    continue
                if node.support < threshold_value:
                    #print node
                    node.delete()
        #print(tree.write())
            out_fd.write(tree.write(format=args.input_tree_format))
            out_fd.write("\n")
            tree_index += 1
