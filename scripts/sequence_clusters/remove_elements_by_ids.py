#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import SequenceClusterRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--cluster_file", action="store", dest="cluster_file", required=True,
                    help="File with clusters")
parser.add_argument("-d", "--black_list_file", action="store", dest="black_list_file", required=True,
                    help="File with ids of elements")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="File to write extracted clusters")
parser.add_argument("-m", "--mode", action="store", dest="mode", default="full",
                    help="Element id comparison mode. Allowed: partial, full(comparison)")

args = parser.parse_args()


SequenceClusterRoutines.remove_elements_by_ids_from_files(args.cluster_file, args.output, args.black_list_file,
                                                          mode=args.mode)
