#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse

import igraph as ig


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input file with graph")
parser.add_argument("-f", "--graph_format", action="store", dest="graph_format",  default='ncol',
                    help="Format of graph. Allowed formats: ncol, lgl, adjacency, dimacs, edgelist, graphviz, gml,"
                         " graphml, graphmlz, pajek, pickle Default: ncol")
parser.add_argument("-l", "--layout", action="store", dest="layout", default='fr',
                    help="Layout used to draw graph")
parser.add_argument("-o", "--output_file", action="store", dest="output", required=True,
                    help="Output file with picture of graph")
args = parser.parse_args()

graph = ig.Graph.Read(args.input, format=args.graph_format)
#graph.es["label"] = map(int, graph.es["weight"])
graph.es["label"] = graph.es["weight"]
graph.vs["label"] = graph.vs["name"]
ig.plot(graph, args.output, layout=args.layout)
