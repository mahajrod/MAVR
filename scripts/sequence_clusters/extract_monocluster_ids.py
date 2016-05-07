#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse

from Routines import SequenceClusterRoutines, FileRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_cluster_dir", action="store", dest="input_cluster_dir", required=True,
                    type=FileRoutines.check_path,
                    help="Directory with files with clusters")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    type=FileRoutines.check_path,
                    help="File to write ids of monoclusters")

args = parser.parse_args()

SequenceClusterRoutines.extract_monocluster_ids_from_file(args.input_cluster_dir, args.output)
