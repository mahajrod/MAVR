#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import sys
import argparse
from RouToolPa.Routines import SequenceClusterRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input syn file")

parser.add_argument("-k", "--key_column_index", action="store", dest="key_column_index", type=int, default=0,
                    help="Index of key column in synonym file. Default: 0")
parser.add_argument("-v", "--value_column_index", action="store", dest="value_column_index", type=int, default=1,
                    help="Index of value column in synonym file.Default: 1")
parser.add_argument("-a", "--allow_repeats_of_key", action="store_true", dest="allow_repeats_of_key", default=False,
                    help="Allow repeats of keys. Default: False")

parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output fam file")


args = parser.parse_args()

SequenceClusterRoutines.syn2fam(args.input, args.output,
                                key_column=args.key_column_index,
                                value_column=args.value_column_index, separator="\t")
