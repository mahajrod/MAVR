#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import SequenceClusterRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="File with clusters")
parser.add_argument("-f", "--label_file", action="store", dest="label_file",
                    help="File with element labels to extract")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="File to write extracted clusters")
parser.add_argument("-p", "--label position", action="store", dest="label_position", default="first",
                    help="Position of label. Allowed - first, last. Default - first")
parser.add_argument("-s", "--separator", action="store", dest="separator", default="@",
                    help="Separator to use. default - '@'")
parser.add_argument("-l", "--label_list", action="store", dest="label_list", type=lambda s: s.split(","),
                    help="Comma-separated list of element labels to extract")

args = parser.parse_args()

if args.label_file and args.label_list:
    raise ValueError("Both --label_file and --label_list were set")
elif (not args.label_file) and (not args.label_list):
    raise ValueError("Neither --label_file or --label_list was set")

SequenceClusterRoutines.extract_clusters_and_elements_by_labels_from_files(args.input,
                                                                           args.label_file if args.label_file else args.label_list,
                                                                           args.output,
                                                                           separator=args.separator,
                                                                           label_position=args.label_position)
