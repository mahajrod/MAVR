#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Parsers.Primer3 import CollectionPrimer3

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with Primer3 output")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-g", "--min_gap_len", action="store", dest="min_gap_len", default=5, type=int,
                    help="Minimum length of polyN to be treated as gap. Default: 5")

args = parser.parse_args()

filtered_results_file = "%s.filtered.res" % args.output_prefix
filtered_out_results_file = "%s.filtered_out.res" % args.output_prefix
primer3_results = CollectionPrimer3(primer3_file=args.input, from_file=True)

primer3_results.remove_primers_with_gaps_in_pcr_product(args.min_gap_len)
primer3_filtered_results, primer3_filtered_out_results = primer3_results.filter_out_records_without_primers()

primer3_filtered_results.write(filtered_results_file)
primer3_filtered_out_results.write(filtered_out_results_file)
