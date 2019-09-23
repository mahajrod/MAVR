#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.Routines import AlignmentRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-c", "--coverage_file_list", action="store", dest="coverage_file_list",
                    type=AlignmentRoutines.split_string_by_comma, required=True,
                    help="Comma-separated list of coverage files")
parser.add_argument("-m", "--mean_coverage_list", action="store", dest="mean_coverage_list",
                    type=lambda s: list(map(int, AlignmentRoutines.split_string_by_comma(s))), required=True,
                    help="Comma-separated list of mean coverage.")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file with mask")
parser.add_argument("-s", "--sample_label_list", action="store", dest="sample_label_list",
                    type=AlignmentRoutines.split_string_by_comma,
                    help="Comma-separated list of sample labels for header in coverage file. Default: use filenames")
parser.add_argument("-x", "--max_coverage_threshold", action="store", dest="max_coverage_threshold", type=float,
                    help="Maximum coverage threshold to treat position as unmasked. Default: not set")
parser.add_argument("-n", "--min_coverage_threshold", action="store", dest="min_coverage_threshold", type=float,
                    help="Minimum coverage threshold to treat position as unmasked. Default: not set")
parser.add_argument("-y", "--min_sample_number", action="store", dest="min_sample_number", type=int,
                    default=1,
                    help="Minimum sample number to mask position. Default: 1")
parser.add_argument("--pos_col", action="store", dest="pos_col", type=int,
                    default=1,
                    help="Index of position column in coverage file (0-based). Default: 1")
parser.add_argument("--cov_col", action="store", dest="cov_col", type=int,
                    default=2,
                    help="Index of coverage column in coverage file (0-based). Default: 2")
parser.add_argument("--scaf_col", action="store", dest="scaf_col", type=int,
                    default=0,
                    help="Index of scaffold column in coverage file (0-based). Default: 0")
args = parser.parse_args()

AlignmentRoutines.calculate_masking_from_coverage_files(args.coverage_file_list, args.mean_coverage_list,
                                                        args.output,
                                                        sample_labels=args.sample_label_list,
                                                        max_threshold=args.max_coverage_threshold,
                                                        min_threshold=args.min_coverage_threshold,
                                                        min_sample_number=args.min_sample_number,
                                                        scaffold_column=args.scaf_col,
                                                        position_column=args.pos_col,
                                                        coverage_column=args.cov_col)
