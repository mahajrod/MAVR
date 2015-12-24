#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

from ete2 import Tree, TreeStyle, AttrFace, faces


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
parser.add_argument("-n", "--feature_name", action="store", dest="feature_name", required=True,
                    help="Name of feature to show")

args = parser.parse_args()


def layout(node):
    #if node.up is not None:

    #attr = AttrFace("decrease", fsize=7, fgcolor="red", text_prefix="-")
    #faces.add_face_to_node(attr, node, 0, position="branch-bottom")
    if node.is_leaf():
        attr = AttrFace("name", fsize=10, fgcolor="black", text_prefix="  ", fstyle="italic")
        faces.add_face_to_node(attr, node, 0, position="aligned")
        attr = AttrFace(args.feature_name, fsize=7, fgcolor="green", text_prefix="")
        faces.add_face_to_node(attr, node, 0, position="branch-top")

#ts.show_branch_support = True
tree_index = 1
with open(args.input_tree_file, "r") as in_fd:
    for line in in_fd:
        tree_line = line.strip()
        tree = Tree(tree_line, format=args.input_tree_format)
        for node in tree.traverse():
            node.img_style["size"] = 0

        ts = TreeStyle()
        ts.layout_fn = layout
        ts.optimal_scale_level = "full"
        ts.branch_vertical_margin = 10
        ts.show_leaf_name = False
        #print(tree.write(format=args.input_tree_format, features=[args.feature_name]))
        if args.draw_mode == "render":
            tree.render("%s_%i.png" % (args.output_tree_file_prefix, tree_index), tree_style=ts,
                        w=args.width, units=args.width_units)
        else:
            tree.show()

        tree_index += 1
