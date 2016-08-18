#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse

from collections import OrderedDict

from Bio import SeqIO

from Routines import SequenceRoutines
from CustomCollections.GeneralCollections import TwoLvlDict

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file_list", action="store", dest="input_file_list", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of files with different assemblies")
parser.add_argument("-l", "--labels_list", action="store", dest="labels_list",
                    type=lambda s: s.split(","),
                    help="Comma-separated list of assembly labels. Should have same length as list of "
                         "input files with assemblies. Default - not set, assemblies will be named like A1, A2, ../ ")
parser.add_argument("-t", "--thresholds", action="store", dest="thresholds", default=[0, 100, 250, 500, 1000],
                    type=lambda s: map(int, s.split(",")),
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
for i in range(0, len(args.input_file_list)):
    assembly_label = args.labels_list[i] if args.labels_list else "A%i" % (i + 1)
    tmp_index = "%s.tmp.idx" % assembly_label
    assemblies_dict[assembly_label] = SeqIO.index_db(tmp_index, args.input_file_list[i],
                                                     format=args.format)

assembly_N50_dict = TwoLvlDict()
assembly_number_of_contigs = TwoLvlDict()
assembly_bins = []
assembly_contig_cumulative_length = TwoLvlDict()
assembly_contig_number_values = TwoLvlDict()

for assembly in assemblies_dict:
    N50_dict, number_of_contigs_dict, total_length, longest_contig, bins, contig_cumulative_length_values, \
        contig_number_values = SequenceRoutines.calculate_assembly_stats(assemblies_dict[assembly],
                                                              thresholds_list=args.thresholds)
    assembly_N50_dict[assembly] = N50_dict
    assembly_number_of_contigs[assembly] = assembly_number_of_contigs
    assembly_contig_cumulative_length[assembly] = contig_cumulative_length_values
    assembly_contig_number_values[assembly] = contig_number_values

    if len(assembly_bins) < len(bins):
        assembly_bins = bins
number_of_bins = len(bins) - 1

# add zeroes to absent bins for all assemblies
for assembly in assembly_contig_cumulative_length:
    bin_number_difference = number_of_bins - len(assembly_contig_cumulative_length[assembly])
    if bin_number_difference > 0:
        assembly_contig_cumulative_length[assembly] += [0 for i in range(0, bin_number_difference)]
        assembly_contig_number_values[assembly] += [0 for i in range(0, bin_number_difference)]

assembly_N50_dict.write("%s.N50" % args.output_prefix)
assembly_number_of_contigs.write("%s.contig_number" % args.output_prefix)
#assembly_bins.write("%s.bins" % args.output_prefix)
assembly_contig_cumulative_length.write("%s.cumulative_length" % args.output_prefix)
assembly_contig_number_values.write("%s.contig_number_values" % args.output_prefix)

for assembly_label in assemblies_dict:
    os.remove("%s.tmp.idx" % assembly_label)


