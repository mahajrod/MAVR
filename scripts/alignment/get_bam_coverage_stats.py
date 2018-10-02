#!/usr/bin/env python
__author__ = 'mahajrod'

import argparse

from Tools.Bedtools import GenomeCov

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_bam", action="store", dest="input_bam", required=True,
                    help="Input bam.")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")


args = parser.parse_args()

GenomeCov.threads = args.threads
GenomeCov.get_bam_coverage_stats(args.input_bam, args.output_prefix, genome_bed=None)

