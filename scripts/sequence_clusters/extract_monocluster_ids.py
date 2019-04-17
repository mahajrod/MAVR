#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import SequenceClusterRoutines, FileRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_cluster_dir", action="store", dest="input_cluster_dir", required=True,
                    type=FileRoutines.check_path,
                    help="Directory with files with clusters")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="File to write ids of monoclusters")
parser.add_argument("-w", "--white_list_ids", action="store", dest="white_list_ids",
                    help="File with ids from white list. ")

args = parser.parse_args()

SequenceClusterRoutines.extract_monocluster_ids_from_file(args.input_cluster_dir, args.output,
                                                          file_with_white_list_ids=args.white_list_ids)
