#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Tools.RepeatMasking import TRF
from RouToolPa.Routines.File import split_filename


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input_file", required=True,
                    help="Input file with sequences")
parser.add_argument("-p", "--path_to_trf", action="store", dest="path_to_trf", default=["", "trf", None],
                    type=split_filename,
                    help="Path to TRF")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", default="trf_report",
                    help="Prefix of output files")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int,
                    help="Number of threads")

args = parser.parse_args()

TRF.threads = args.threads
TRF.path = args.path_to_trf[0]
TRF.cmd = args.path_to_trf[1] + (args.path_to_trf[2] if args.path_to_trf[2] else "")

TRF.parallel_search_long_tandem_repeat(args.input_file, args.output_prefix,
                                       report_flanking_sequences=True, splited_fasta_dir="splited_fasta_dir",
                                       splited_result_dir="splited_output", converted_output_dir="converted_output",
                                       max_len_per_file=100000, store_intermediate_files=False, max_repeat_length=6)
