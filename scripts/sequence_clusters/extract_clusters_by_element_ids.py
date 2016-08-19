#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse

from Routines import SequenceClusterRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--cluster_file", action="store", dest="cluster_file", required=True,
                    help="File with clusters")
parser.add_argument("-d", "--element_file", action="store", dest="element_file", required=True,
                    help="File with ids of elements")
parser.add_argument("-o", "--output", action="store", dest="output", default="stdout",
                    help="File to write extracted clusters")
parser.add_argument("-m", "--mode", action="store", dest="mode", default="w",
                    help="extraction mode. Allowed - 'w' - if elements from element_id_list are present "
                         "in cluster extracts only that elements 'a' - if elements from element_id_list are present "
                         "in cluster extracts all elements Default - 'w'")
args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

SequenceClusterRoutines.extract_clusters_by_element_ids_from_file(args.cluster_file, args.element_file, out_fd,
                                                                  mode=args.mode)
