#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
import os

from Tools.RepeatMasking import TRF

from Routines.File import split_filename


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input_file", required=True,
                    help="Input file with sequences")
parser.add_argument("-p", "--path_to_trf", action="store", dest="path_to_trf", default=["", "trf", None],
                    type=split_filename,
                    help="Path to TRF")
parser.add_argument("-o", "--output_file", action="store", dest="output_file", default="trf_report.gff",
                    help="Output gff file")
parser.add_argument("-m", "--matching_weight", action="store", dest="matching_weight", default=2, type=int,
                    help="Matching weight. Default: 2")
parser.add_argument("-s", "--mismatching_penalty", action="store", dest="mismatching_penalty", default=7, type=int,
                    help="Mismatching penalty. Default: 7")
parser.add_argument("-l", "--indel_penalty", action="store", dest="indel_penalty", default=7, type=int,
                    help="Indel_penalty. Default: 7")
parser.add_argument("-a", "--matching_probability", action="store", dest="matching_probability", default=80, type=int,
                    help="Matching probability. Default: 80")
parser.add_argument("-d", "--indel_probability", action="store", dest="indel_probability", default=10, type=int,
                    help="Indel probability. Default: 10")
parser.add_argument("-c", "--min_score", action="store", dest="min_score", default=50, type=int,
                    help="Minimum alignment score to report. Default: 50")
parser.add_argument("-e", "--max_period_size", action="store", dest="max_period_size", default=500, type=int,
                    help="Maximum period size to report. Default: 500")
parser.add_argument("-f", "--report_flanking_sequence", action="store_true", dest="report_flanking_sequence",
                    help="Report flanking sequence. Default: False")
args = parser.parse_args()

TRF.path = args.path_to_trf[0]
TRF.cmd = args.path_to_trf[1] + (args.path_to_trf[2] if args.path_to_trf[2] else "")

TRF.search_tandem_repeats(args.input_file, matching_weight=args.matching_weight,
                          mismatching_penalty=args.mismatching_penalty, indel_penalty=args.indel_penalty,
                          match_probability=args.matching_probability, indel_probability=args.indel_probability,
                          min_alignment_score=args.min_score, max_period=args.max_period_size,
                          report_flanking_sequences=args.report_flanking_sequences, make_dat_file=True)

trf_report = "%s.%i.%i.%i.%i.%i.%i.%i.dat" % (args.input_file, args.matching_weight, args.mismatching_penalty,
                                              args.indel_penalty, args.matching_probability, args.indel_probability,
                                              args.min_score, args.max_period_size)
TRF.convert_trf_report_to_gff(trf_report, args.output_file)



