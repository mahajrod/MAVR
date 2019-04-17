#!/usr/bin/env python
__author__ = 'mahajrod'

import argparse
from RouToolPa.Tools.Annotation import Exonerate
from RouToolPa.Routines.File import split_filename


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input_file", required=True,
                    help="Input file with sequences")
parser.add_argument("-d", "--database_file", action="store", dest="database_file", required=True,
                    help="File with database")
parser.add_argument("-p", "--path_to_exonerate", action="store", dest="path_to_exonerate", default=["", "exonerate", None],
                    type=split_filename, help="Path to exonerate")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int,
                    help="Number of threads")
parser.add_argument("-m", "--model", action="store", dest="model",
                    help="Model to use")
parser.add_argument("-n", "--number_of_sequences_per_file", action="store", dest="num_of_sequences_per_file",
                    type=int, default=1000,
                    help="Number of sequences per file")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", default="trf_report",
                    help="Prefix of output files")
parser.add_argument("-r", "--store_intermediate_file", action="store_true", dest="store",
                    help="Dont remove intermediate files")
parser.add_argument("--results_to_report", action="store", dest="num_of_results", type=int,
                    help="Number of results to report per query")

args = parser.parse_args()

Exonerate.threads = args.threads
Exonerate.path = args.path_to_exonerate[0]
Exonerate.cmd = args.path_to_exonerate[1] + (args.path_to_exonerate[2] if args.path_to_exonerate[2] else "")

Exonerate.parallel_alignment(args.input_file, args.database_file, args.model,
                             num_of_recs_per_file=args.num_of_sequences_per_file,
                             show_alignment=True, show_sugar=True, show_cigar=None,
                             show_vulgar=None, show_query_gff=None, show_target_gff=True,
                             store_intermediate_files=True,
                             splited_fasta_dir="splited_fasta_dir", splited_result_dir="splited_output",
                             converted_output_dir="converted_output",
                             number_of_results_to_report=args.num_of_results)




