#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines.File import make_list_of_path_to_files
from RouToolPa.Routines import SequenceRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file_list", action="store", dest="input", required=True,
                    type=lambda s: make_list_of_path_to_files(s.split(",")),
                    help="Comma-separated list of input files/directories with sequences")
parser.add_argument("-o", "--output_file", action="store", dest="output", default="stdout",
                    help="Output file with renamed sequences")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of input and output files. Allowed formats genbank, fasta(default)")
parser.add_argument("-s", "--syn_file", action="store", dest="syn_file", required=True,
                    help="File to write correspondences between new and old ids")
parser.add_argument("-n", "--numerical_part_of_id_length", action="store", dest="numerical_part_of_id_length",
                    default=8, type=int,
                    help="Length of numerical part of id. Default: 8")
parser.add_argument("-d", "--id_prefix", action="store", dest="id_prefix", required=True,
                    help="Prefix of new sequence ids")
parser.add_argument("-l", "--clear_description", action="store_true", dest="clear_description", default=False,
                    help="Clear description. Default - False")

args = parser.parse_args()

SequenceRoutines.rename_records_by_sequential_ids_from_files(args.input, args.output, args.syn_file, format=args.format,
                                                             clear_description=args.clear_description,
                                                             record_id_prefix=args.id_prefix,
                                                             length_of_numerical_part=args.numerical_part_of_id_length,
                                                             parse_mode="parse", index_file="temp.idx")

