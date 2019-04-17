#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import FastQRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input fastq file")
parser.add_argument("-o", "--output_dir", action="store", dest="output_dir", required=True,
                    help="Output directory")
parser.add_argument("-p", "--output_prefix", action="store", dest="output_prefix",
                    help="Prefix of output files. Default: not set")
parser.add_argument("-s", "--output_suffix", action="store", dest="output_suffix", default=".fastq",
                    help="Suffix of output files. Default: .fastq")

args = parser.parse_args()

FastQRoutines.split_illumina_fastq_by_lanes(args.input, args.output_dir,
                                            output_prefix=args.output_prefix,
                                            output_suffix=args.output_suffix)
