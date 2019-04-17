#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import FastQRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input fastq file")
parser.add_argument("-o", "--output_file", action="store", dest="output", required=True,
                    help="output fastq file")

args = parser.parse_args()

FastQRoutines.reverse_complement(args.input, args.output)
