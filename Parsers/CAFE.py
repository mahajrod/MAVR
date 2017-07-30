#!/usr/bin/env python
__author__ = 'mahajrod'

import numpy as np
from ete2 import Tree, TreeStyle, AttrFace, faces

from Parsers.Abstract import Header


class HeaderCAFE(list, Header):
    def __str__(self):
        return "#" + "\t".join(self)


class MetadataCAFE():
    def __init__(self, format1, format2):
        self.format1 = format1
        self.format2 = format2


class RecordCAFE():
    def __init__(self, id, tree, family_p_value):
        self.id = id
        self.tree = tree
        self.family_p_value = family_p_value

    def __str__(self):
        return "%s\t%f\t%s" % (self.id, self.family_p_value,
                               self.tree.write(format=100,
                                               features=self.tree.features - set(["support", "name", "dist"]),
                                               format_root_node=True))

    def node_str(self):
        node_string = "%s\t%f" % (self.id, self.family_p_value)
        features = self.tree.features - set(["support", "name", "dist"])
        for node in self.tree.traverse():
            list_of_features = [eval("node.%s" % feature) for feature in features]
            list_of_features = map(lambda x: str(x) if x is not None else ".", list_of_features)
            node_string += "\t" + ",".join(list_of_features)
        return node_string


class GeneralDataCAFE():
    def __init__(self, tree, lambda_value=None):
        self.tree = tree
        for node in self.tree.traverse():
            node.img_style["size"] = 0
        self.lambda_value = lambda_value
        self.node_id_list = [node.id for node in self.tree.traverse()]

    def __str__(self):
        tree_str = "##Common tree\t%s\n" % self.tree.write(format=8,
                                                           features=self.tree.features - set(["support", "name"]),
                                                           format_root_node=True)
        lambda_str = "##Lambda\t%f" % self.lambda_value
        return tree_str + lambda_str

    def draw(self, out_file_prefix, w=140, units="mm"):
        def layout(node, feature):
            if node.up is not None or feature == "id":
                attr = AttrFace(feature, fsize=5, fgcolor="green")
                faces.add_face_to_node(attr, node, 0, position="branch-top")
        for feature in self.tree.features - set(["support", "name"]):
            if feature[0] == "_":
                continue            # skip _nid feature and other technical features of tree

            def layout_arg(node):
                return layout(node, feature)
            ts = TreeStyle()
            ts.layout_fn = layout_arg
            self.tree.render("%s_%s.png" % (out_file_prefix, feature), w=w, units=units, tree_style=ts)

    def draw_expansion_contraction(self, outfile_prefix="expansion_contraction_tree"):

        tree = self.tree
        for node in tree.traverse():
            node.img_style["size"] = 0

        def layout(node):
            if node.up is not None:
                attr = AttrFace("expansion", fsize=7, fgcolor="green", text_prefix="+")
                faces.add_face_to_node(attr, node, 0, position="branch-top")
                attr = AttrFace("decrease", fsize=7, fgcolor="red", text_prefix="-")
                faces.add_face_to_node(attr, node, 0, position="branch-bottom")
            if node.is_leaf():
                attr = AttrFace("name", fsize=10, fgcolor="black", text_prefix="  ", fstyle="italic")
                faces.add_face_to_node(attr, node, 0, position="aligned")
        ts = TreeStyle()
        ts.layout_fn = layout
        ts.optimal_scale_level = "full"
        ts.branch_vertical_margin = 10
        ts.show_leaf_name = False
        #ts.allow_face_overlap =True
        tree.render("%s.png" % outfile_prefix, w=200, units='mm', tree_style=ts, dpi=1200)
        tree.render("%s.svg" % outfile_prefix, w=200, units='mm', tree_style=ts, dpi=1200)

    def draw_significant_expansion_contraction(self, outfile_prefix="significant_expansion_contraction_tree"):

        tree = self.tree
        for node in tree.traverse():
            node.img_style["size"] = 0

        def layout(node):
            if node.up is not None:
                attr = AttrFace("significant_expansion", fsize=7, fgcolor="green", text_prefix="+")
                faces.add_face_to_node(attr, node, 0, position="branch-top")

                attr = AttrFace("significant_contraction", fsize=7, fgcolor="red", text_prefix="-")
                faces.add_face_to_node(attr, node, 0, position="branch-bottom")

            if node.is_leaf():
                attr = AttrFace("name", fsize=10, fgcolor="black", text_prefix="  ", fstyle="italic")
                faces.add_face_to_node(attr, node, 0, position="aligned")
        ts = TreeStyle()
        ts.layout_fn = layout
        ts.optimal_scale_level = "full"
        ts.branch_vertical_margin = 10
        ts.show_leaf_name = False
        #ts.allow_face_overlap =True
        tree.render("%s.png" % outfile_prefix, w=200, units='mm', tree_style=ts, dpi=1200)
        tree.render("%s.svg" % outfile_prefix, w=200, units='mm', tree_style=ts, dpi=1200)

    def write_general_tree(self, out_file):
        with open(out_file, "w") as out_fd:
            out_fd.write(self.tree.write(format=8,
                                         features=self.tree.features - set(["support", "name"]),
                                         format_root_node=True))
            out_fd.write("\n")


class ReportCAFE():
    def __init__(self, report_file=None, from_file=True, records=None, metadata=None, general_data=None, header=None):
        if from_file:
            with open(report_file, "r") as report_fd:
                common_tree = Tree(report_fd.readline().strip()[5:] + ";")
                lambda_value = float(report_fd.readline().strip().split(":")[-1])
                id_tree = Tree(report_fd.readline().strip().split(":")[-1] + ";", format=8)
                for id_tree_node, common_tree_node in zip(id_tree.traverse(), common_tree.traverse()):
                    node_id = int(id_tree_node.name.split("<")[-1][:-1])
                    common_tree_node.add_feature("id", node_id)
                data_format1 = self.split_format1(report_fd.readline().strip().split(":")[-1])
                data_format2 = report_fd.readline().strip().split(":")[-1].split()
                data_format2[0] = data_format2[0][1:-1]
                for i in range(1, len(data_format2)):
                    data_format2[i] = data_format2[i][:-1]
                data_format2 = map(lambda x: int(x), data_format2)
                for data_type in [float, int, int, int]:
                    name, data = report_fd.readline().strip().split(":")
                    data = self.split_format1(data, data_type=data_type)
                    name = name.strip().replace(" ", "_").lower()
                    for node in common_tree.traverse():
                        node.add_feature(name, data[data_format1.index(node.id)] if node.id in data_format1 else None)

                self.general_data = GeneralDataCAFE(common_tree, lambda_value=lambda_value)
                self.metadata = MetadataCAFE(data_format1, data_format2)
                self.header = HeaderCAFE(report_fd.readline().strip().split("\t"))
                self.records = []
                for line in report_fd:
                    id, tree, p_value, vitebri_p_value = line.strip().split()
                    p_value = float(p_value)
                    #print(tree)
                    tree = Tree(tree + ";", format=1) # was 8
                    #print(vitebri_p_value)
                    vitebri_p_value = self.split_format1(vitebri_p_value[1:-1], separator=",", data_type=float)
                    for tree_node, id_node in zip(tree.traverse(), self.general_data.tree.traverse()):
                        if tree_node.is_leaf:
                            tree_node.name, number_of_genes = tree_node.name.split("_")
                        else:
                            number_of_genes = tree_node.name[1:]
                        #print (node.name)
                        number_of_genes = int(number_of_genes)
                        tree_node.add_feature("N", number_of_genes)
                        tree_node.add_feature("p_value", vitebri_p_value[data_format1.index(id_node.id)] if id_node.id in data_format1 else None)
                        parent_tree_node = tree_node.up
                        if parent_tree_node is None:
                            delta = None
                        else:
                            delta = tree_node.N - parent_tree_node.N
                        tree_node.add_feature("delta", delta)

                    self.records.append(RecordCAFE(id, tree, p_value))
        else:
            self.general_data = general_data
            self.metadata = metadata
            self.header = header
            self.records = records

    def __iter__(self):
        for record in self.records:
            yield record

    def split_format1(self, format1_string, separator=None, data_type=int):

        tmp_format1 = format1_string.replace(")", "").replace("(", "").split(separator)
        format1 = []
        #print tmp_format1
        for i in range(0, len(tmp_format1)):
            format1 += map(lambda x: None if x == '-' else data_type(x), tmp_format1[i].split(","))
        return format1

    def convert(self, out_file, output_type="tree"):
        # ouput_type - possible values: tree, node
        with open(out_file, "w") as out_fd:
            out_fd.write(str(self.general_data))
            out_fd.write("\n")
            if output_type == "tree":
                header = "#id\tfamily p-value\ttree\n"
                out_fd.write(header)
                for record in self:
                    out_fd.write(str(record))
                    out_fd.write("\n")
            elif output_type == "node":
                out_fd.write("#Node_format: %s\n" % (",".join(self.records[0].tree.features - set(["support", "name", "dist"]))))
                header = "#id\tfamily_p_value"
                for node in self.general_data.tree.traverse():
                    header += "\t%s" % node.id
                header += "\n"
                out_fd.write(header)
                for record in self:
                    out_fd.write(record.node_str())
                    out_fd.write("\n")

    def get_per_node_values(self, filter_by_p_value=False, p_value_cutoff=0.05):
        node_values_dict = dict((node_id, []) for node_id in self.general_data.node_id_list)
        feature_names_list = ["family_id"] + [feature for feature in self.records[0].tree.features - set(["support", "name", "dist"])]
        for record in self:
            node_index = 0
            for node in record.tree.traverse():
                if not filter_by_p_value or node.p_value <= p_value_cutoff:
                    node_values_dict[self.general_data.node_id_list[node_index]].append([record.id] + [eval("node.%s" % feature) for feature in node.features - set(["support", "name", "dist"])])
                node_index += 1
        return node_values_dict, feature_names_list

    #----------------------------General filtering functions--------------------------------------
    def filter_records(self, expression_function):
        #expression_function must return  boolean value
        filtered = []
        filtered_out = []
        for record in self:
            if expression_function(record):
                filtered.append(record)
            else:
                filtered_out.append(record)
        return filtered, filtered_out

    def filter(self, expression_function):
        filtered_records, filtered_out_records = self.filter_records(expression_function)

        return ReportCAFE(from_file=False, records=filtered_records, general_data=self.general_data,
                          metadata=self.metadata, header=self.header), \
               ReportCAFE(from_file=False, records=filtered_out_records, general_data=self.general_data,
                          metadata=self.metadata, header=self.header),
    #---------------------------------------------------------------------------------------------

    #--------------------------Custom filtering functions-----------------------------------------
    def filter_by_family_p_value(self, p_value_cutoff):
        return self.filter(lambda record: record.family_p_value <= p_value_cutoff)


if __name__ == "__main__":
    cafe_report_file = "/home/mahajrod/genetics/Dobrzhansky/project/gene_families/data/pashas/resultfilereport.cafe"
    cafe_report = ReportCAFE(cafe_report_file)
    #print(cafe_report.general_data.tree.write(features=["id", "remain", "Expansion", "decrease", "average_expansion"]))
    #print(cafe_report.general_data.lambda_value)
    #print(cafe_report.records[3].tree.write(features=["N", "delta", "p_value"], format=8, format_root_node=True))

    cafe_report.general_data.tree.render("mytree.png", w=183, units="mm")
    cafe_report.general_data.draw("labeled_tree")
