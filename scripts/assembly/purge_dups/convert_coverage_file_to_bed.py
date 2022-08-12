#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse

from RouToolPa.Tools.AsseblyQC import PurgeDups

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with purge_dups coverage estimations. Usually it is 'TX.base.cov'")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
args = parser.parse_args()

PurgeDups.convert_coverage_file_to_bed(args.input_file, args.output_prefix)

