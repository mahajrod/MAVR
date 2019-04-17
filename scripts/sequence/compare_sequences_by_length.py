#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import SequenceRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-a", "--seq_file_a", action="store", dest="seq_file_a", required=True,
                    help="Sequence file A")
parser.add_argument("-b", "--seq_file_b", action="store", dest="seq_file_b", required=True,
                    help="Sequence file B")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")


args = parser.parse_args()


SequenceRoutines.compare_sequences_by_length_from_file(args.seq_file_a, args.seq_file_b, args.output_prefix)


