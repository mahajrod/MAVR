#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import SequenceClusterRoutines, FileRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with clusters")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="File to write selected_clusters")
parser.add_argument("-w", "--white_list_ids", action="store", dest="white_list_ids",
                    help="File with ids from white list. ")
parser.add_argument("-n", "--min", action="store", dest="min", type=int,
                    help="Minimal number of elements in cluster. Default: not set")
parser.add_argument("-x", "--max", action="store", dest="max", type=int,
                    help="Maximal number of elements in cluster. Default: not set")
args = parser.parse_args()

SequenceClusterRoutines.extract_clusters_by_size_from_file(args.input,
                                                           min_cluster_size=args.min,
                                                           max_cluster_size=args.max,
                                                           white_list_ids=args.white_list_ids,
                                                           out_file=args.output)
