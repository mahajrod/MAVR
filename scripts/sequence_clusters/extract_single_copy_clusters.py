#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import SequenceClusterRoutines, FileRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Comma-separated list of files/directories with clusters")
parser.add_argument("-l", "--label", action="store_true", dest="label",
                    help="Label elements by file basenames. Default - False")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="File to write clusters with single-copy clusters")
parser.add_argument("-p", "--label position", action="store", dest="label_position", default="first",
                    help="Position of label. Allowed - first, last. Default - first")
parser.add_argument("-s", "--separator", action="store", dest="separator", default="@",
                    help="Separator to use. default - '@'")

args = parser.parse_args()

list_of_cluster_files = FileRoutines.make_list_of_path_to_files(args.input)

single_copy_clusters = SequenceClusterRoutines.extract_single_copy_clusters_from_files(list_of_cluster_files, args.output,
                                                                label_elements=args.label, separator=args.separator,
                                                                label_position=args.label_position)

print("Was found %i single-copy clusters" % len(single_copy_clusters))