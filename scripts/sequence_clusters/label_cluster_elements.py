#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse

from Routines import SequenceClusterRoutines, FileRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--cluster_file", action="store", dest="cluster_file", required=True,
                    help="File with clusters")
parser.add_argument("-l", "--label", action="store", dest="label",
                    help="Label to use. If not set a basename of file is used")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="File to write clusters with labeled elements")
parser.add_argument("-p", "--label position", action="store", dest="label_position", default="first",
                    help="Position of label. Allowed - first, last. Default - first")
parser.add_argument("-s", "--separator", action="store", dest="separator", default="@",
                    help="Separator to use. default - '@'")

args = parser.parse_args()

label = args.label if args.label else FileRoutines.split_filename(args.cluster_file)[1]

SequenceClusterRoutines.label_cluster_elements_from_file(args.cluster_file, label, args.output,
                                                         separator=args.separator,
                                                         label_position=args.label_position)
