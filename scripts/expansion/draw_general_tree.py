#!/usr/bin/env python
import argparse
from ete2 import Tree, TreeStyle, AttrFace, faces


def layout(node):
    if node.up is not None:
        attr = AttrFace("expansion", fsize=7, fgcolor="green", text_prefix="+")
        faces.add_face_to_node(attr, node, 0, position="branch-top")
        attr = AttrFace("decrease", fsize=7, fgcolor="red", text_prefix="-")
        faces.add_face_to_node(attr, node, 0, position="branch-bottom")
    if node.is_leaf():
        attr = AttrFace("name", fsize=10, fgcolor="black", text_prefix="  ", fstyle="italic")
        faces.add_face_to_node(attr, node, 0, position="aligned")


def layout_significant(node):
    if node.up is not None:
        attr = AttrFace("significant_expansion", fsize=7, fgcolor="green", text_prefix="+")
        faces.add_face_to_node(attr, node, 0, position="branch-top")
        attr = AttrFace("significant_contraction", fsize=7, fgcolor="red", text_prefix="-")
        faces.add_face_to_node(attr, node, 0, position="branch-bottom")
    if node.is_leaf():
        attr = AttrFace("name", fsize=10, fgcolor="black", text_prefix="  ", fstyle="italic")
        faces.add_face_to_node(attr, node, 0, position="aligned")


def layout_simple(node):
    if node.is_leaf():
        attr = AttrFace("name", fsize=10, fgcolor="black", text_prefix="  ", fstyle="italic")
        faces.add_face_to_node(attr, node, 0, position="aligned")

parser = argparse.ArgumentParser()

parser.add_argument("-c", "--tree_with_general_data", action="store", dest="tree_with_general_data",
                    help="Tree with general data about expansion analysis")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    help="Prefix of output file with drawn trees")

args = parser.parse_args()

common_tree_file = args.tree_with_general_data

general_tree_file = "%s.png" % args.output_prefix
general_significant_tree_file = "%s_significant.png" % args.output_prefix
general_tree_file_no_labels = "%s_no_labels.png" % args.output_prefix
with open(common_tree_file, "r") as tree_fd:
    tree = Tree(tree_fd.readline().strip(), format=1)
        
for node in tree.traverse():
    node.img_style["size"] = 0

ts = TreeStyle()
ts.layout_fn = layout
ts.optimal_scale_level = "full"
ts.branch_vertical_margin = 10
ts.show_leaf_name = False
#ts.allow_face_overlap =True

tree.render(general_tree_file, w=200, units='mm', tree_style=ts, dpi=300)

ts_significant = TreeStyle()
ts_significant.layout_fn = layout_significant
ts_significant.optimal_scale_level = "full"
ts_significant.branch_vertical_margin = 10
ts_significant.show_leaf_name = False
#ts.allow_face_overlap =True

tree.render(general_significant_tree_file, w=200, units='mm', tree_style=ts_significant, dpi=300)

ts_simple = TreeStyle()
ts_simple.layout_fn = layout_simple
#ts_simple.optimal_scale_level = "full"
ts_simple.branch_vertical_margin = 15
ts_simple.show_leaf_name = False


tree.render(general_tree_file_no_labels, w=200, units='mm', tree_style=ts_simple, dpi=300)