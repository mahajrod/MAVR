#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import FileRoutines, SequenceClusterRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--cluster_id_file", action="store", dest="cluster_id_file",
                    help="File with ids of clusters to extract. Extract all clusters if not set")
parser.add_argument("-f", "--cluster_file", action="store", dest="cluster_file", required=True,
                    help="File with clusters")
parser.add_argument("-p", "--seq_file", action="store", dest="seq_file", required=True,
                    type=lambda s: FileRoutines.make_list_of_path_to_files(s.split(",")),
                    help="List of comma-separated files/directories with sequences")
parser.add_argument("-r", "--seq_file_format", action="store", dest="seq_file_format", default="fasta",
                    help="Format of file with sequences")
parser.add_argument("-c", "--create_dir_for_each_cluster", action="store_true", dest="create_dir_for_each_cluster",
                    help="Create separate directory for each cluster")
parser.add_argument("-o", "--output_prefix", action="store", dest="output",
                    help="Output prefix to use")
parser.add_argument("-d", "--output_directory", action="store", dest="out_dir", default="./",
                    help="Directory to write output")
parser.add_argument("-n", "--dont_skip_cluster_if_absent_element", action="store_true",
                    dest="dont_skip_cluster_if_absent_element", default=False,
                    help="Don't skip cluster with absent sequences")


args = parser.parse_args()

SequenceClusterRoutines.extract_sequences_from_selected_clusters(args.cluster_id_file, args.cluster_file,
                                                                 args.seq_file, output_dir=args.out_dir,
                                                                 seq_format=args.seq_file_format,
                                                                 out_prefix=args.output,
                                                                 create_dir_for_each_cluster=args.create_dir_for_each_cluster,
                                                                 skip_cluster_if_no_sequence_for_element=not args.dont_skip_cluster_if_absent_element,
                                                                 parsing_mode="parse")
