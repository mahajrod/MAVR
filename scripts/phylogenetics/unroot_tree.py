#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import TreeRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_tree_file", action="store", dest="input_tree_file", required=True,
                    help="File with input trees to unroot")
parser.add_argument("-o", "--output_tree_file", action="store", dest="output_tree_file", required=True,
                    help="File with output unrooted trees")
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
100 	topology only.
Default: 0""")
parser.add_argument("-t", "--output_tree_format", action="store", dest="output_tree_format", type=int, default=None,
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
100 	topology only.
Default: same as input format""")

args = parser.parse_args()

TreeRoutines.unroot_tree_from_file(args.input_tree_file, args.output_tree_file,
                                   input_tree_format=args.input_tree_format, output_tree_format=args.output_tree_format)
