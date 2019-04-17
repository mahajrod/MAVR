#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import SequenceClusterRoutines, FileRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="File with input clusters")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="File to write output_clusters")
parser.add_argument("-s", "--syn_file", action="store", dest="syn_file",
                    help="File with synonyms for labels. If not set only changing position of label is possible")

parser.add_argument("-p", "--input_label position", action="store", dest="input_label_position", default="first",
                    help="Position of label in input file. Allowed - first, last. Default - first")
parser.add_argument("-s", "--input_separator", action="store", dest="input_separator", default="@",
                    help="Separator between label and element id in input file. Default - '@'")

parser.add_argument("-a", "--output_label position", action="store", dest="output_label_position", default="first",
                    help="Position of label in output file. Allowed - first, last. Default - first")
parser.add_argument("-r", "--output_separator", action="store", dest="output_separator", default="@",
                    help="Separator between label and element id in output file. Default - '@'")

args = parser.parse_args()

label = args.label if args.label else FileRoutines.split_filename(args.cluster_file)[1]

SequenceClusterRoutines.replace_label_from_file(args.input_file, args.output_file, args.syn_file,
                                                old_separator=args.input_separator,
                                                old_label_position=args.input_label_position,
                                                new_separator=args.output_separator,
                                                new_label_position=args.output_label_position)
