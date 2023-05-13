#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import MGIRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-s", "--sample_file", action="store", dest="sample_file", required=True,
                    help="Sample file")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Output prefix")
parser.add_argument("-m", "--max_index_len", action="store", dest="max_index_len", default=10,
                    help="maximum allowed length of index. Default: 10")
parser.add_argument("-e", "--max_index_errors", action="store", dest="max_index_errors", default=2, type=int,
                    help="Maximum number of errors allowed per each index. Default: 2")

args = parser.parse_args()

sample_df = MGIRoutines.parse_samples(args.sample_file)
MGIRoutines.analyze_sample_df(sample_df, args.output_prefix, max_index_len=args.max_index_len,
                              max_index_errors=args.max_index_errors)
