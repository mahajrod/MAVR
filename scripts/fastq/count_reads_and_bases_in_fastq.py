#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import FastQRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    type=FastQRoutines.make_list_of_path_to_files_from_string,
                    help="Comma-separated list of fastq files")
parser.add_argument("-o", "--output", action="store", dest="output",
                    help="Output file. Default: stdout only")

args = parser.parse_args()

FastQRoutines.count_reads_and_bases(args.input, args.output)
