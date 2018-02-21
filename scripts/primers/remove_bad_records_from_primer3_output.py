#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse

from Parsers.Primer3 import CollectionPrimer3

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with Primer3 output")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file.")
parser.add_argument("-g", "--min_gap_len", action="store", dest="min_gap_len", default=5, type=int,
                    help="Minimum length of polyN to be treated as gap. Default: 5")

args = parser.parse_args()

primer3_results = CollectionPrimer3(primer3_file=args.input, from_file=True)

primer3_results.remove_primers_with_gaps_in_pcr_product(args.min_gap_len)
primer3_results.filter_out_records_without_primers()
primer3_results.write(args.output)
