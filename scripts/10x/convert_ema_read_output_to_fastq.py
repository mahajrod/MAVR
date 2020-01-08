#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.Tools.Alignment import EMA


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input file with reads in EMA format. Default: stdin")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

args = parser.parse_args()

EMA.convert_ema_read_format_to_fastq(args.input, args.output_prefix, bufferig=1000000000)
