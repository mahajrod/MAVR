#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Tools.RepeatMasking import TRF



parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input gff with TRF repeats")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output gff with filtered repeats")
parser.add_argument("-x", "--filtered_out", action="store", dest="filtered_out", required=True,
                    help="Gff with filtered out repeats")
parser.add_argument("-n", "--min_period", action="store", dest="min_period", type=int,
                    help="Minimum period of repeats to extract")
parser.add_argument("-m", "--max_period", action="store", dest="max_period", type=int,
                    help="Maximum period of repeats to extract")
parser.add_argument("-b", "--min_copy_number", action="store", dest="min_copy_number", type=float,
                    help="Minimum number of copies(float) to extract")
parser.add_argument("-c", "--max_copy_number", action="store", dest="max_copy_number", type=float,
                    help="Maximum number of copies(float) to extract")

parser.add_argument("-p", "--pattern", action="store", dest="pattern",
                    help="Extract patterns only with this pattern")
parser.add_argument("-d", "--min_percentage_of_matches", action="store", dest="min_percentage_of_matches", type=int,
                    help="Minimum percentage of matches(int)")
parser.add_argument("-e", "--max_percentage_of_indels", action="store", dest="max_percentage_of_indels", type=int,
                    help="Maximum percentage of indels(int)")

"""
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int,
                    help="Number of threads")
parser.add_argument("-p", "--path_to_trf", action="store", dest="path_to_trf", default=["", "trf", None],
                    type=split_filename,
                    help="Path to TRF")
"""

args = parser.parse_args()

#TRF.threads = args.threads
#TRF.path = args.path_to_trf[0]
#TRF.cmd = args.path_to_trf[1] + (args.path_to_trf[2] if args.path_to_trf[2] else "")

TRF.filter_trf_gff(args.input, args.output, args.filtered_out, min_period=args.min_period, max_period=args.max_period,
                   min_copy_number=args.min_copy_number, max_copy_number=args.max_copy_number,
                   pattern=args.pattern, min_percentage_of_matches=args.min_percentage_of_matches,
                   max_percentage_of_indels=args.max_percentage_of_indels, min_entropy=None, max_entropy=None)

