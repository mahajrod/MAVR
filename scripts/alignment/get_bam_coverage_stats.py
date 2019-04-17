#!/usr/bin/env python
__author__ = 'mahajrod'

import argparse
from RouToolPa.Tools.Bedtools import GenomeCov

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_bam", action="store", dest="input_bam", required=True,
                    help="Input bam.")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-n", "--min_coverage", action="store", dest="min_coverage", type=float,
                    help="Minimum coverage to use. Default - not set")
parser.add_argument("-m", "--max_coverage", action="store", dest="max_coverage", type=float,
                    help="Maximum coverage to use. Default - not set")

args = parser.parse_args()

#GenomeCov.threads = args.threads
GenomeCov.get_bam_coverage_stats(args.input_bam, args.output_prefix, genome_bed=None,
                                 max_coverage=args.max_coverage, min_coverage=args.min_coverage,
                                 verbose=True)

