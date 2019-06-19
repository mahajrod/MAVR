#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import SequenceRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input fasta")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output fasta")
parser.add_argument("-s", "--syn_file", action="store", dest="syn_file", required=True,
                    help="File with id synonyms")
parser.add_argument("-k", "--key", action="store", dest="key", default=0, type=int,
                    help="Key column(0-based) in synonym file. Default: 0")
parser.add_argument("-v", "--value", action="store", dest="value", default=1, type=int,
                    help="Value column(0-based) in synonym file. Default: 1")

args = parser.parse_args()

SequenceRoutines.make_fasta_by_correspondence(args.input, args.syn_file, args.output,
                                              key_column=args.key, value_column=args.value)
