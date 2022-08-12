#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse

from RouToolPa.Tools.AsseblyQC import PurgeDups

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with purge_dups coverage estimations. Usually it is 'TX.base.cov'")
parser.add_argument("-l", "--length_file", action="store", dest="length_file", required=True,
                    help="File with scaffold_lengths")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file")
args = parser.parse_args()


PurgeDups.add_lengths_to_dups_bed(args.input, args.length_file, args.output)
