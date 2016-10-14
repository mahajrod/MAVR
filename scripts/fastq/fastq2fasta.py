#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

from BConverters import SequenceConverters

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input fastq file")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output fastq file ")

args = parser.parse_args()

SequenceConverters.fastq2fasta(args.input, args.output)
