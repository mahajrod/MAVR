#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

from ete3 import Tree, TreeStyle


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_tree_file", action="store", dest="input_tree_file",
                    help="File with input trees")
parser.add_argument("-o", "--output_tree_file_prefix", action="store", dest="output_tree_file_prefix",
                    help="Prefix of file with output trees")
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

parser.add_argument("-w", "--width", action="store", dest="width", type=int,
                    help="Tree image width")

parser.add_argument("-u", "--width_units", action="store", dest="width_units", default="mm",
                    help="Tree image width units: 'px': pixels, 'mm': millimeters, 'in': inches. Default: 'mm' ")
parser.add_argument("-m", "--draw_mode", action="store", dest="draw_mode", default="render",
                    help="Draw mode - render to file or show preview. Default - render")

args = parser.parse_args()


#ts.show_branch_support = True
tree_index = 1
with open(args.input_tree_file, "r") as in_fd:
    for line in in_fd:
        tree_line = line.strip()
        tree = Tree(tree_line, format=args.input_tree_format)
        for node in tree.traverse():
            node.img_style["size"] = 0
        if args.draw_mode == "render":
            tree.render("%s_%i.png" % (args.output_tree_file_prefix, tree_index), w=args.width, units=args.width_units)
        else:
            tree.show()

        tree_index += 1