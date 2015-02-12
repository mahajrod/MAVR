#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from collections import OrderedDict

from ete2 import Tree


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_tree_file", action="store", dest="input_tree_file",
                    help="File with input trees")
parser.add_argument("-o", "--output_tree_file_prefix", action="store", dest="output_tree_file_prefix",
                    help="Prefix of file with otput trees")
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
parser.add_argument("-m", "--min_children_tree_size", action="store", dest="min_children_tree_size", type=int,
                    default=1, help="Minimum size of children tree. Default: 1")
parser.add_argument("-s", "--stat_file_prefix", action="store", dest="stat_file_prefix", default="statistics",
                    help="Prefix of file with size (number of leaves in tree) of children trees")

args = parser.parse_args()

tree_index = 1


with open(args.input_tree_file, "r") as in_fd:
    for line in in_fd:
        children_index = 1
        stat_dict = OrderedDict({})
        tree_line = line.strip()
        tree = Tree(tree_line, format=args.input_tree_format)
        total = 0
        for children_tree in tree.children:
            if len(children_tree) < args.min_children_tree_size:
                continue
            children_tree_name = "%s_%i_%i" % (args.output_tree_file_prefix, tree_index, children_index)
            stat_dict[children_tree_name] = len(children_tree)
            total += stat_dict[children_tree_name]
            with open("%s.nwk" % children_tree_name, "w") as out_fd:
                out_fd.write(children_tree.write(format=args.input_tree_format))
            children_index += 1

        with open("%s_tree_%i.t" % (args.stat_file_prefix, tree_index), "w") as stat_fd:
            stat_fd.write("#Totaly_trees\t%i\n" % len(stat_dict))
            stat_fd.write("#Totaly_leaves\t%i\n" % total)
            stat_fd.write("#Tree\tNumber_of_leaves\n")
            for children_tree in stat_dict:
                stat_fd.write("%s\t%i\n" % (children_tree, stat_dict[children_tree]))
        tree_index += 1