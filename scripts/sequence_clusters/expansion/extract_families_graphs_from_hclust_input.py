#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os
import argparse
from multiprocessing import Pool
from RouToolPa.Collections.General import SynDict
from RouToolPa.Routines.File import check_path


"""

def call_hcluster(input, output, edge_weight, edge_density):
    hcluster_options = " -w %i -s %f" % (edge_weight, edge_density)
    hcluster_string = "%s %s %s > %s" % (args.hcluster_sg_path, hcluster_options, input, output)
    print (hcluster_string)
    os.system(hcluster_string)


def call_hcluster_single_arg(arg_list):
    call_hcluster(args.input, *arg_list)
"""

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--hclust_input_file", action="store", dest="hclust_input", required=True,
                    help="File with input graph for hclust")
parser.add_argument("-c", "--families_file", action="store", dest="fam_file", required=True,
                    help="File with families")
parser.add_argument("-o", "--output_dir", action="store", dest="output_dir",
                    help="Directory to write output")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads to use")
args = parser.parse_args()

args.output_dir = check_path(args.output_dir)


def check_edge_strict(nodes_list, id_list):
    for node in nodes_list:
        if node not in id_list:
            return False
    return True


def check_edge_soft(nodes_list, id_list):
    for node in nodes_list:
        if node in id_list:
            return True
    return False

families_dict = SynDict()
families_dict.read(args.fam_file, separator="\t", split_values=True, values_separator=",")

try:
    os.mkdir(args.output_dir)
except OSError:
    pass

graph_list = []
with open(args.hclust_input, "r") as in_fd:
    for line in in_fd:
        graph_list.append(line.strip().split("\t"))

def extract_fam_graph(family_name):
    print("Started extraction for family %s" % family_name)
    family_genes_ids = families_dict[family_name]
    try:
        os.mkdir("%s%s" % (args.output_dir, family_name))
    except OSError:
        pass

    fam_soft_fd = open("%s%s/%s_with_outer_edges.graph" % (args.output_dir, family_name, family_name), "w")
    """
    with open(args.hclust_input, "r") as in_fd:
        for line in in_fd:
            edge_nodes = line.split("\t")[:2]
            if check_edge_soft(edge_nodes, family_genes_ids):
                fam_soft_fd.write(line)
    """
    for edge in graph_list:
        if check_edge_soft(edge[:-1], family_genes_ids):
            fam_soft_fd.write("\t".join(edge) + "\n")
    fam_soft_fd.close()
    fam_strict_fd = open("%s%s/%s.graph" % (args.output_dir, family_name, family_name), "w")
    with open("%s%s/%s_with_outer_edges.graph" % (args.output_dir, family_name, family_name), "r") as in_fd:
        for line in in_fd:
            edge_nodes = line.split("\t")[:2]
            if check_edge_strict(edge_nodes, family_genes_ids):
                fam_strict_fd.write(line)
    fam_strict_fd.close()

pool = Pool(args.threads)
pool.map(extract_fam_graph, families_dict.keys())

