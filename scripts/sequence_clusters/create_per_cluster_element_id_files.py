#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import SequenceClusterRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_cluster_file", action="store", dest="input_cluster_file", required=True,
                    help="Input file with clusters")

parser.add_argument("-c", "--cluster_column_index", action="store", dest="cluster_column_index", type=int, default=0,
                    help="Index of cluster column in synonym file. Default: 0")
parser.add_argument("-v", "--element_column_index", action="store", dest="element_column_index", type=int, default=1,
                    help="Index of element column in synonym file.Default: 1")
parser.add_argument("-e", "--separator", action="store", dest="column_separator", default='\t',
                    help="Column separator in synonym file. Default: \\t")

parser.add_argument("-o", "--output_directory", action="store", dest="output_dir", required=True,
                    help="Output directory")

args = parser.parse_args()

SequenceClusterRoutines.create_per_cluster_element_id_files_from_file(args.input_cluster_file,
                                                                      args.output_dir,
                                                                      cluster_column=args.cluster_column_index,
                                                                      element_column=args.element_column_index,
                                                                      column_separator=args.column_separator,
                                                                      element_separator=",")
