#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Collections.General import IdList


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--hclust_input_file", action="store", dest="hclust_input", required=True,
                    help="File with input graph for hclust")
parser.add_argument("-d", "--id_file", action="store", dest="id_file", required=True,
                    help="File with ids of genes(nodes) to remove from graph")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    help="Prefix of output files")
parser.add_argument("-m", "--mode", action="store", dest="mode", default="one",
                    help="Removal mode. Possible variants: both, one. Default: one")
args = parser.parse_args()

id_list = IdList()
id_list = id_list.read(args.id_file)

filtered_fd = open("%s_filtered.t" % args.output_prefix, "w")
filtered_out_fd = open("%s_filtered_out.t" % args.output_prefix, "w")

if args.mode == "one":
    def expression(nodes_list):
        for node in nodes_list:
            if node in id_list:
                return False
        return True
elif args.mode == "both":
    def expression(nodes_list):
        for node in nodes_list:
            if node not in id_list:
                return True
        return False

with open(args.hclust_input, "r") as in_fd:
    for line in in_fd:
        edge_nodes = line.split("\t")[:2]
        if expression(edge_nodes):
            filtered_fd.write(line)
        else:
            filtered_out_fd.write(line)

filtered_out_fd.close()
filtered_fd.close()