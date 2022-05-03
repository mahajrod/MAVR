#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

from RouToolPa.Routines import FastQRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-f", "--forward_reads", action="store", dest="forward_reads", required=True,
                    help="Input file with forward reads")
parser.add_argument("-r", "--reverse_reads", action="store", dest="reverse_reads", required=True,
                    help="Input file with reverse reads")
parser.add_argument("-x", "--readname_prefix", action="store", dest="readname_prefix", required=True,
                    help="Prefix to use for read names")
parser.add_argument("-n", "--interleaved", action="store_true", dest="interleaved", default=False,
                    help="Output renamed reads in interleaved format")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files ")

args = parser.parse_args()

FastQRoutines.add_read_name_prefix(args.forward_reads, args.reverse_reads, args.readname_prefix,
                                   args.output_prefix, interleaved=args.interleaved)
