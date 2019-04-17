#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import ExpressionRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-c", "--count_table_file", action="store", dest="count_table_file", required=True,
                    help="File with count table")
parser.add_argument("-l", "--transcript_length_file", action="store", dest="transcript_length_file",
                    help="File with transcript lengths")
parser.add_argument("-o", "--output_file", action="store", dest="output", required=True,
                    help="Output file with FPKM counts")

args = parser.parse_args()

ExpressionRoutines.calculate_fpkm_for_count_table(args. count_table_file, args.transcript_length_file, args.output,
                                                  separator="\t")
