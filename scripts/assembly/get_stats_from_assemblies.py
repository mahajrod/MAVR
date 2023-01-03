#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
import multiprocessing as mp
from collections import OrderedDict

import pandas as pd

from RouToolPa.Parsers.Sequence import CollectionSequence


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file_list", action="store", dest="input_file_list", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of files with different assemblies")
parser.add_argument("-l", "--labels_list", action="store", dest="labels_list",
                    type=lambda s: s.split(","),
                    help="Comma-separated list of assembly labels. Should have same length as list of "
                         "input files with assemblies. Default - not set, assemblies will be named like A1, A2, ../ ")
parser.add_argument("-e", "--thresholds", action="store", dest="thresholds", default=[0, 100, 250, 500, 1000],
                    type=lambda s: list(map(int, s.split(","))),
                    help="Comma-separated list of thresholds for N50 calculations. "
                         "Default: 0,100,250,500,1000")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of input files")

args = parser.parse_args()

if args.labels_list is not None:
    if len(args.labels_list) != len(args.input_file_list):
        raise ValueError("Length of labels list is not equal to number of files with assemblies")

assemblies_dict = OrderedDict()
stats_dict = OrderedDict({"N50": pd.DataFrame(),
                          "L50": pd.DataFrame(),
                          "Ns": pd.DataFrame(),
                          "Total length": pd.DataFrame(),
                          "Total scaffolds": pd.DataFrame()})
thresholds_stats = OrderedDict([(threshold, pd.DataFrame()) for threshold in args.thresholds])

for i in range(0, len(args.input_file_list)):
    assembly_label = args.labels_list[i] if args.labels_list else "A%i" % (i + 1)
    assembly_output_prefix = "%s.%s" % (args.output_prefix, assembly_label)
    assemblies_dict[assembly_label] = CollectionSequence(in_file=args.input_file_list[i], parsing_mode="parse").get_stats_and_features(thresholds_list=args.thresholds, count_gaps=True)
    assemblies_dict[assembly_label].to_csv("%s.tsv" % assembly_output_prefix, sep="\t")
    print(assemblies_dict[assembly_label])
    for stat_entry in stats_dict:
        stats_dict[stat_entry][assembly_label] = assemblies_dict[assembly_label].loc[stat_entry]

    for threshold in thresholds_stats:
        thresholds_stats[threshold][assembly_label] = assemblies_dict[assembly_label][threshold]

for stat_entry in stats_dict:
    stats_dict[stat_entry].to_csv("%s.%s.tsv" % (args.output_prefix, stat_entry.replace(" ", "_")), sep="\t")

for threshold in thresholds_stats:
    thresholds_stats[threshold].to_csv("%s.%i.tsv" % (args.output_prefix, threshold), sep="\t")
