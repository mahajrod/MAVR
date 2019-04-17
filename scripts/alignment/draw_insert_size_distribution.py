#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Tools.Samtools import SamtoolsV1


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input sam/bam file")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    help="Prefix of output files")
parser.add_argument("-s", "--separator", action="store", dest="separator", default="\n",
                    help="Separator between values in input file. Default - '\\n', i.e. one value per line")

parser.add_argument("-w", "--width_of_bins", action="store", dest="width_of_bins", type=float, default=5,
                    help="Width of bins in histogram. Default: 5")
parser.add_argument("-n", "--min_value", action="store", dest="min_insert_size", type=float, default=0,
                    help="Minimum insert to show on histogram. Default - 0")
parser.add_argument("-x", "--max_value", action="store", dest="max_insert_size", type=float, default=1000,
                    help="Maximum value to show. Default: 1000")
parser.add_argument("-g", "--logbase", action="store", dest="logbase", type=int, default=10,
                    help="Logbase to use for log-scaled histograms")
parser.add_argument("-e", "--extensions", action="store", dest="extensions", type=lambda x: x.split(","),
                    default=["png"],
                    help="Comma-separated list of extensions for histogram files. Default: png only")

args = parser.parse_args()

SamtoolsV1.draw_insert_size_distribution(args.input, args.output_prefix, width_of_bin=args.width_of_bins,
                                         max_insert_size=args.max_insert_size,
                                         min_insert_size=args.min_insert_size, extensions=args.extensions,
                                         separator=args.separator, logbase=args.logbase)
