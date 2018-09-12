#!/usr/bin/env python
__author__ = 'mahajrod'
import os
import argparse

from Routines import AnnotationsRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-c", "--chunk_dir", action="store", dest="chunk_dir", required=True,
                    help="Input dir with chunks")
parser.add_argument("-s", "--chunk_suffix", action="store", dest="chunk_suffix", default=".out",
                    help="Suffix of chunk files. Chunks must be named as <prefix>_<chunk number><suffix> . "
                         "Default: '.out'")
parser.add_argument("-n", "--total_number_of_chunks", action="store", dest="number_of_chunks", type=int,
                    required=True,
                    help="Total number of chunks")
parser.add_argument("-m", "--min_chunk_size", action="store", dest="min_chunk_size", type=int,
                    required=True,
                    help="Minimum size of chunk file.")

args = parser.parse_args()

AnnotationsRoutines.check_chunks(args.chunk_dir, args.number_of_chunks, args.min_chunk_size,
                                 chunk_name_suffix=args.chunk_suffix)
