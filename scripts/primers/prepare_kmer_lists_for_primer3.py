#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Kmers import Glistmaker
from RouToolPa.Routines import FileRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated_list_of_input_files")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads. Default: 1")

parser.add_argument("-p", "--glistmaker_path", action="store", dest="glistmaker_path",
                    type=FileRoutines.split_filename, default=["", "glistmaker", ""],
                    help="Path to Glistmaker binary")

args = parser.parse_args()

Glistmaker.path = args.glistmaker_path[0]
Glistmaker.cmd = args.glistmaker_path[1] + args.glistmaker_path[2]
Glistmaker.threads = args.threads

Glistmaker.generate_kmer_lists_for_primer3(args.input, args.output_prefix)
