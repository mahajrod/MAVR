#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from collections import OrderedDict
from ete2 import Tree
from RouToolPa.Collections.General import SynDict

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_tree_file", action="store", dest="input_tree_file",  required=True,
                    help="File with input trees to collapse")
parser.add_argument("-o", "--output_tree_file", action="store", dest="output_tree_file",
                    help="File to write tree with features")
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

parser.add_argument("-e", "--feature_file", action="store", dest="feature_file", required=True,
                    help="File with features")
parser.add_argument("-s", "--species_column", action="store", dest="species_column", type=int, default=0,
                    help="Number of column(0-based) with species names. Default: 1")
parser.add_argument("-a", "--feature_columns", action="store", dest="feature_columns",
                    type=lambda s: map(int, s.split(",")), default=[1],
                    help="Comma-separated list of numbers of columns(0-based) with "
                         "features in feature file. Default: 1")
parser.add_argument("-n", "--feature_names", action="store", dest="feature_names",
                    type=lambda s: s.split(","), default=["feature"],
                    help="Comma-separated list of names of features from feature file "
                         "Number of feature names should be the same as number of columns "
                         "set by -a/--feature_columns option. Default: feature")
parser.add_argument("-c", "--comments_prefix", action="store", dest="comments_prefix",
                    help="Prefix of comments in feature file. Line with this prefix will be ignored")

args = parser.parse_args()

if len(args.feature_names) != len(args.feature_columns):
    raise ValueError("Number of feature names is not equal to number of feature columns")
features_dict = OrderedDict()

for feature_name, feature_column in zip(args.feature_names, args.feature_columns):
    features_dict[feature_name] = SynDict()
    # print args.species_column,  feature_column
    features_dict[feature_name].read(args.feature_file, header=False, separator="\t",
                                     split_values=False, values_separator=",", key_index=args.species_column,
                                     value_index=feature_column,
                                     comments_prefix=args.comments_prefix)

tree_index = 1
with open(args.input_tree_file, "r") as in_fd:
    with open(args.output_tree_file, "w") as out_fd:
        for line in in_fd:
            tree_line = line.strip()
            tree = Tree(tree_line, format=args.input_tree_format)
            for node in tree.traverse():
                for feature in features_dict:
                    if node.name in features_dict[feature]:
                        node.add_feature(feature, features_dict[feature][node.name])
            out_fd.write(tree.write(format=args.input_tree_format, features=args.feature_names))
            out_fd.write("\n")
            tree_index += 1
