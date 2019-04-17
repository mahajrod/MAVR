#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import sys
import argparse
from RouToolPa.Routines import FastQRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input fastq file")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output file with tiles. Default: stdout")

args = parser.parse_args()

FastQRoutines.find_tiles(args.input, args.output)
