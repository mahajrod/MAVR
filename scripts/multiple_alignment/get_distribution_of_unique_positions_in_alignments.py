#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from collections import OrderedDict
from RouToolPa.Routines import MultipleAlignmentRoutines, FileRoutines, MatplotlibRoutines
from RouToolPa.Collections.General import TwoLvlDict


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    type=lambda s: FileRoutines.make_list_of_path_to_files(s.split(",")),
                    help="Comma-separated list of files/directories with alignments")
parser.add_argument("-o", "--output_directory", action="store", dest="output_dir", default="./",
                    help="Output directory to write count files. Default - current directory")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of alignments")
parser.add_argument("-g", "--gap_symbol", action="store", dest="gap_symbol", default="-",
                    help="Gap symbol. Default - '-'")

parser.add_argument("-m", "--histogram_output", action="store", dest="histogram_output", required=True,
                    help="File to write histogram")

args = parser.parse_args()

unique_position_dict = TwoLvlDict()

FileRoutines.safe_mkdir(args.output_dir)

for alignment_file in args.input:
    alignment_name_list = FileRoutines.split_filename(alignment_file)
    output_prefix = "%s/%s.unique_positions" % (args.output_dir, alignment_name_list[1])

    unique_position_dict[alignment_name_list[1]] = MultipleAlignmentRoutines.count_unique_positions_per_sequence_from_file(alignment_file,
                                                                                                                           output_prefix,
                                                                                                                           format=args.format,
                                                                                                                           gap_symbol="-",
                                                                                                                           return_mode="relative",
                                                                                                                           verbose=False)

species_list = unique_position_dict.sl_keys()

data_dict = OrderedDict()

for species in species_list:
    data_dict[species] = []
    for alignment in unique_position_dict:
        data_dict[species].append(unique_position_dict[alignment][species])

data_list = [data_dict[species] for species in data_dict]

MatplotlibRoutines.extended_percent_histogram(data_list, args.histogram_output, input_mode="percent", label=species_list)
