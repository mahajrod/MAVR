#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import SequenceClusterRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_cluster_file", action="store", dest="input_cluster_file", required=True,
                    help="Input file with clusters")
parser.add_argument("-s", "--syn_file", action="store", dest="syn_file", required=True,
                    help="File with synonyms to elements")

parser.add_argument("-k", "--key_column_index", action="store", dest="key_column_index", type=int, default=0,
                    help="Index of key column in synonym file. Default: 0")
parser.add_argument("-v", "--value_column_index", action="store", dest="value_column_index", type=int, default=1,
                    help="Index of value column in synonym file.Default: 1")
parser.add_argument("-e", "--separator", action="store", dest="column_separator", default='\t',
                    help="Column separator in synonym file. Default: \\t")

parser.add_argument("-o", "--output_cluster_file", action="store", dest="output_cluster_file", required=True,
                    help="File to write clusters with renamed elements")
parser.add_argument("-a", "--elements_without_synonyms_file", action="store", dest="elements_without_synonyms_file",
                    help="File to write cluster elements without synonyms. Default: don't write")
parser.add_argument("-r", "--remove_clusters_with_not_renamed_elements", action="store_true",
                    dest="remove_clusters_with_not_renamed_elements",
                    help="Remove clusters with not renamed elements. Default: false ")

args = parser.parse_args()

SequenceClusterRoutines.rename_elements_in_clusters(args.input_cluster_file, args.syn_file, args.output_cluster_file,
                                                    remove_clusters_with_not_renamed_elements=args.remove_clusters_with_not_renamed_elements,
                                                    syn_file_key_column_index=args.key_column_index,
                                                    syn_file_value_column_index=args.value_column_index,
                                                    syn_file_column_separator=args.column_separator,
                                                    elements_with_absent_synonyms_file=args.elements_without_synonyms_file)
