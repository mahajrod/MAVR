#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import SequenceClusterRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--cluster_file", action="store", dest="cluster_file", required=True,
                    help="File with clusters")
parser.add_argument("-d", "--element_file", action="store", dest="element_file", required=True,
                    help="File with ids of elements")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="File to write extracted clusters")
parser.add_argument("-c", "--id_column_index", action="store", dest="id_column_index", type=int,
                    help="Index(0-based) of id column in id file. ")
parser.add_argument("-m", "--mode", action="store", dest="mode", default="w",
                    help="extraction mode. Allowed - 'w' - if elements from element_id_list are present "
                         "in cluster extracts only that elements; 'a' - if elements from element_id_list are present "
                         "in cluster extracts all elements. Default - 'w'")

args = parser.parse_args()


SequenceClusterRoutines.extract_clusters_by_element_ids_from_file(args.cluster_file, args.element_file, args.output,
                                                                  mode=args.mode, id_column=args.id_column_index)
