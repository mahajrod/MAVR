#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse

from RouToolPa.Routines import AlignmentRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-c", "--coverage_bed", action="store", dest="coverage_bed",
                    required=True,
                    help="Bed file with coverage")
parser.add_argument("-m", "--mean_coverage", action="store", dest="mean_coverage", type=float,
                    help="Mean coverage to use for mask calculation.")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file with mask")
parser.add_argument("--cov_col", action="store", dest="cov_col", type=int,
                    default=3,
                    help="Index of coverage column in coverage file (0-based). Default: 3")

# relative threshold flags
parser.add_argument("-x", "--max_coverage_threshold", action="store", dest="max_coverage_threshold", type=float,
                    help="Maximum coverage threshold to treat position as unmasked. Default: not set")
parser.add_argument("-n", "--min_coverage_threshold", action="store", dest="min_coverage_threshold", type=float,
                    help="Minimum coverage threshold to treat position as unmasked. Default: not set")

# absolute threshold flags
parser.add_argument("--max_coverage_absolute", action="store", dest="max_coverage_absolute",
                    type=float, help="Absolute max coverage threshold to treat position as unmasked. Default: not set")
parser.add_argument("--min_coverage_absolute", action="store", dest="min_coverage_absolute",
                    type=float, help="Absolute min coverage threshold to treat position as unmasked. Default: not set")

args = parser.parse_args()

# Convert absolute -> relative if needed
if args.max_coverage_absolute:
    max_threshold = args.max_coverage_absolute / args.mean_coverage
else:
    max_threshold = args.max_coverage_threshold

if args.min_coverage_absolute:
    min_threshold = args.min_coverage_absolute / args.mean_coverage
else:
    min_threshold = args.min_coverage_threshold

AlignmentRoutines.calculate_masking_from_coverage_bed(
    args.coverage_bed,
    args.mean_coverage,
    args.output,
    max_threshold=max_threshold,
    min_threshold=min_threshold,
    coverage_column=args.cov_col,
    buffering=100000000
)
