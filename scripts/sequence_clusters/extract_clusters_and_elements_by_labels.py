#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse

from Routines import SequenceClusterRoutines, FileRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="File with clusters")
parser.add_argument("-l", "--label_file", action="store", dest="label_file", required=True,
                    help="File with labels to extract")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="File to write extracted_clusters")
parser.add_argument("-p", "--label position", action="store", dest="label_position", default="first",
                    help="Position of label. Allowed - first, last. Default - first")
parser.add_argument("-s", "--separator", action="store", dest="separator", default="@",
                    help="Separator to use. default - '@'")

args = parser.parse_args()

SequenceClusterRoutines.extract_clusters_and_elements_by_labels_from_files(args.input, args.label_file, args.output,
                                                                           separator=args.separator,
                                                                           label_position=args.label_position)
