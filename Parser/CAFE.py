#!/usr/bin/env python
__author__ = 'mahajrod'

from ete2 import Tree

from Parser.Abstract import Header


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

class GeneralDataCAFE():
    def __init__(self, tree, lambda_value=None):
        self.tree = tree
        self.lambda_value = lambda_value


class ReportCAFE():
    def __init__(self, report_file):
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
                tree = Tree(tree + ";", format=8)
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

    def split_format1(self, format1_string, separator=None, data_type=int):

        tmp_format1 = format1_string.replace(")", "").replace("(", "").split(separator)
        format1 = []
        #print tmp_format1
        for i in range(0, len(tmp_format1)):
            format1 += map(lambda x: None if x == '-' else data_type(x), tmp_format1[i].split(","))
        return format1


if __name__ == "__main__":
    cafe_report_file = "/home/mahajrod/genetics/Dobrzhansky/project/gene_families/data/pashas/resultfilereport.cafe"
    cafe_report = ReportCAFE(cafe_report_file)
    #print(cafe_report.general_data.tree.write(features=["id", "remain", "expansion", "decrease", "average_expansion"]))
    #print(cafe_report.general_data.lambda_value)
    print(cafe_report.records[3].tree.write(features=["N", "delta", "p_value"], format=8))
    for node in cafe_report.records[3].tree.traverse():
        print(node.name, node.N)
        parent_node = node.up
        if parent_node is None:
            print("None")
        else:
            print(parent_node.name, parent_node.N)
        print ("")