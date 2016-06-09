#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from BConverters import SequenceConverters

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input fasta file")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output phylip file")

args = parser.parse_args()

SequenceConverters.fasta2phylip(args.input, args.output)

