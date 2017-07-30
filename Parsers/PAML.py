#!/usr/bin/env python

import sys
from collections import OrderedDict

from ete2 import Tree #, TreeStyle, AttrFace, faces, NodeStyle
from CustomCollections.GeneralCollections import IdList, TwoLvlDict, SynDict


class CodeMLReport():
    """
    Only parsing  of  dN, dS and omega trees is supported
    """
    def __init__(self, reportfile, treefile=None):
        self.dStree = None
        self.dNtree = None
        self.Wtree = None
        self.LnL = None
        with open(reportfile, "r") as in_fd:
            for line in in_fd:
                if line[:8] == "dS tree:":
                    self.dStree = Tree(in_fd.next())
                    continue
                elif line[:8] == "dN tree:":
                    self.dNtree = Tree(in_fd.next())
                    continue
                elif line.strip() == "w ratios as labels for TreeView:":
                    self.Wtree = Tree(in_fd.next().replace(" #", ": "))  # convert omega tree to readable by ete2 format
                    continue
                elif line[:10] == "lnL(ntime:":
                    self.LnL = float(line.strip().split(":")[-1].split()[0])
                    continue
        if treefile:
            with open(treefile, "r") as tree_fd:
                self.tree = Tree(tree_fd.readline())
            i = 1
            if self.dStree and self.dNtree and self.Wtree:
                for tree_node, dNnode, dSnode, Wnode in zip(self.tree.traverse(), self.dNtree.traverse(), self.dStree.traverse(), self.Wtree.traverse()):
                    tree_node.add_feature("dN", dNnode.dist)
                    tree_node.add_feature("dS", dSnode.dist)
                    tree_node.add_feature("W", Wnode.dist)
                    tree_node.add_feature("id", i)
                    i += 1
        else:
            self.tree = None

    def write_trees(self, treefile_prefix, format=5):
        dNtree_file = "%s_dN.nwk" % treefile_prefix
        dStree_file = "%s_dS.nwk" % treefile_prefix
        Wtree_file = "%s_W.nwk" % treefile_prefix
        tree_file = "%s_all.nwk" % treefile_prefix
        with open(dNtree_file, "w") as dn_fd:
            dn_fd.write(self.dNtree.write(format=format) + "\n")
        with open(dStree_file, "w") as ds_fd:
            ds_fd.write(self.dStree.write(format=format) + "\n")
        with open(Wtree_file, "w") as w_fd:
            w_fd.write(self.Wtree.write(format=format) + "\n")
        if self.tree:
            with open(tree_file, "w") as t_fd:
                t_fd.write(self.tree.write(format=format, features={"dN", "dS", "W"}) + "\n")

    def get_all_values(self, output_file):
        with open(output_file, "w") as out_fd:
            out_fd.write("#dN\tdS\tW\n")
            for dNnode, dSnode, Wnode in zip(self.dNtree.traverse(), self.dStree.traverse(), self.Wtree.traverse()):
                if dNnode.is_root():
                    continue
                out_fd.write("%f\t%f\t%f\n" % (dNnode.dist, dSnode.dist, Wnode.dist))

    @staticmethod
    def _get_tree_dist_values(tree, expression=None):
        feature_values_list = IdList()
        for node in tree.traverse():
            if expression is not None:
                if not expression(node):
                    continue
            feature_values_list.append(node.dist)
        return feature_values_list

    @staticmethod
    def _get_tree_dist_dict(tree):
        feature_values_dict = OrderedDict()
        for node in tree.traverse():
            if not node.is_leaf():
                continue
            feature_values_dict[node.name] = node.dist
        return feature_values_dict

    def get_leaf_values(self, write=True):
        leaf_values_dict = TwoLvlDict()
        dN_dict = self._get_tree_dist_dict(self.dNtree)
        dS_dict = self._get_tree_dist_dict(self.dStree)
        W_fict = self._get_tree_dist_dict(self.Wtree)

        leaf_values_dict["dN"] = dN_dict
        leaf_values_dict["dS"] = dS_dict
        leaf_values_dict["W"] = W_fict

        if write:
            leaf_values_dict.write("leaf_values.t")
        return leaf_values_dict

    def find_leaves_with_positive_selection(self, write=True):
        leaf_values_dict = self.get_leaf_values(write=False)
        positive_selected_leaves_dict = SynDict()
        for leaf_name in leaf_values_dict["W"]:
            if leaf_values_dict["W"][leaf_name] > 1:
                positive_selected_leaves_dict[leaf_name] = leaf_values_dict["W"][leaf_name]
        if write:
            positive_selected_leaves_dict.write("leaves_with_positive_selection.t")
        return positive_selected_leaves_dict

    def get_feature_values(self, mode="all"):
        if mode == "all":
            expression = lambda node: not node.is_root()
            suffix = "_all.t"
        elif mode == "leaves":
            expression = lambda node: node.is_leaf() and (not node.is_root())
            suffix = "_leaves.t"
        elif mode == "internal":
            expression = lambda node: not node.is_leaf() and (not node.is_root())
            suffix = "_internal.t"
        else:
            raise ValueError("Possible modes: all, leaves, internal")

        dN_list = self._get_tree_dist_values(self.dNtree, expression=expression)
        dS_list = self._get_tree_dist_values(self.dStree, expression=expression)
        W_list = self._get_tree_dist_values(self.Wtree, expression=expression)

        dN_list.write("dN%s" % suffix)
        dS_list.write("dS%s" % suffix)
        W_list.write("W%s" % suffix)

        return dN_list, dS_list, W_list

    """
    def draw_trees(self, output_prefix, pic_width=900):
        tree = self.tree
        print self.tree
        print tree.write(format=5, features={"dN", "dS", "W"})
        red_style = NodeStyle()
        red_style["fgcolor"] = "#0f0f0f"
        red_style["size"] = 0
        red_style["vt_line_color"] = "#ff0000"
        red_style["hz_line_color"] = "#ff0000"
        red_style["vt_line_width"] = 8
        red_style["hz_line_width"] = 8

        for node in tree.traverse():
            node.img_style["size"] = 0
            node.dS = "%.2f" % node.dS
            node.dN = "%.2f" % node.dN
            node.dW = "%.2f" % node.dW
            #if node.W > 1.0:
            #    node.set_style(red_style)

        def layout(node):
            if node.up is not None:
                attr = AttrFace("dN", fsize=4, fgcolor="green", text_prefix="")
                faces.add_face_to_node(attr, node, 0, position="branch-top")
                attr = AttrFace("dS", fsize=4, fgcolor="red", text_prefix="")
                faces.add_face_to_node(attr, node, 0, position="branch-bottom")
            if node.is_leaf():
                attr = AttrFace("name", fsize=10, fgcolor="black", text_prefix="  ", fstyle="italic")
                faces.add_face_to_node(attr, node, 0, position="aligned")

        def omega_layout(node):
            if node.up is not None:
                attr = AttrFace("W", fsize=4, fgcolor="black", text_prefix="")
                faces.add_face_to_node(attr, node, 0, position="branch-top")

            if node.is_leaf():
                attr = AttrFace("name", fsize=10, fgcolor="black", text_prefix="  ", fstyle="italic")
                faces.add_face_to_node(attr, node, 0, position="aligned")

        ts = TreeStyle()
        ts.layout_fn = layout
        ts.optimal_scale_level = "full"
        ts.branch_vertical_margin = 10
        ts.show_leaf_name = False
        tree.render("%s.png" % output_prefix)
        tree.render("%s_dN_dS_tree.png" % output_prefix, w=pic_width, units='mm', tree_style=ts, dpi=200)

        omega_ts = TreeStyle()
        omega_ts.layout_fn = omega_layout
        omega_ts.optimal_scale_level = "full"
        omega_ts.branch_vertical_margin = 10
        omega_ts.show_leaf_name = False

        tree.render("%s_W_tree.png" % output_prefix, w=pic_width, units='mm', tree_style=omega_ts, dpi=300)
    """
    def branches_with_positive_selection(self):
        for node in self.Wtree.traverse():
            if node.dist > 1.0:
                return True
        return False

    def convert_trees_to_tsv(self, output_prefix):

        w_tab_file = "%s.dn.ds.w.tab" % output_prefix
        id_tree_file = "%s.id.nwk" % output_prefix
        with open(w_tab_file, 'w') as w_fd:
            w_fd.write("#node_id\tnode_name\tdN\tdS\tW\n")
            if self.tree:
                for node in self.tree.traverse():
                    w_fd.write("%i\t%s\t%f\t%f\t%f\n" % (node.id,
                                                         node.name if node.name else ".",
                                                         node.dN,
                                                         node.dS,
                                                         node.W))
                    with open(id_tree_file, "w") as t_fd:
                        t_fd.write(self.tree.write(format=5, features={"dN", "dS", "W", "id"}) + "\n")
            else:
                i = 1
                for dNnode, dSnode, Wnode in zip(self.dNtree.traverse(), self.dStree.traverse(), self.Wtree.traverse()):

                    for node in dNnode, dSnode, Wnode:
                        node.add_feature("id", i)

                    w_fd.write("%i\t%s\t%f\t%f\t%f\n" % (node.id,
                                                         node.name if node.name else ".",
                                                         dNnode.dist,
                                                         dSnode.dist,
                                                         Wnode.dist))
                    i += 1

                    with open(id_tree_file, "w") as t_fd:
                        t_fd.write(self.dNtree.write(format=format, features={"id"}) + "\n")

    def extract_trees(self, out_prefix):
        self.write_trees(out_prefix)
        self.get_feature_values(mode="all")
        self.get_feature_values(mode="leaves")
        self.get_feature_values(mode="internal")
        self.get_all_values("dN_dS_W.t")
        self.get_leaf_values()
        if self.branches_with_positive_selection():
            sys.stderr.write("Presence of branches with positive selection\n")

        self.convert_trees_to_tsv(out_prefix)

