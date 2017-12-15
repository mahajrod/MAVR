#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import argparse

from Tools.Filter import Clumpify

parser = argparse.ArgumentParser()

parser.add_argument("-f", "--forward_reads", action="store", dest="forward_reads", required=True,
                    help="File with forward reads")
parser.add_argument("-r", "--reverse_reads", action="store", dest="reverse_reads", default=None,
                    help="File with reverse reads")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-g", "--gzip_output", action="store_true", dest="gzip_output", default=False,
                    help="gzip output. Default: False")

args = parser.parse_args()

Clumpify.remove_optical_duplicates(args.forward_reads, args.output_prefix,
                                   reverse_reads=args.reverse_reads, gzip_output=args.gzip_output)

