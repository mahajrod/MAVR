#!/usr/bin/env python
__author__ = 'mahajrod'

import sys
import argparse
from pathlib import Path
import pandas as pd
from RouToolPa.Parsers.Sequence import CollectionSequence
from RouToolPa.Routines import SequenceRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input fasta file with sequences. Default: stdin")
parser.add_argument("-d", "--id_list", action="store", dest="id_list",
                    help="Comma-separated list (or file with one id per line) "
                         "of ids of scaffolds to reverse complement."
                         " If not set all scaffolds will be reverse complemented")

parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output fasta file with reverse complement sequences. Default:  stdout")

args = parser.parse_args()

id_series = pd.to_csv(args.id_list, header=None).squeeze("columns") if Path(args.id_list).exists() else pd.Series(args.id_list.split(","))

seq_col = CollectionSequence(in_file=args.input)

for seq_id in id_series:
    seq_col.records[seq_id] = SequenceRoutines.reverse_complement(seq_col.records[seq_id])

seq_col.write(args.output, max_symbols_per_line=60)
