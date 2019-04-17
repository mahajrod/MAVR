#!/usr/bin/env python
__author__ = 'mahajrod'
import argparse
from RouToolPa.Routines import AnnotationsRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-d", "--chunk_dir", action="store", dest="chunk_dir", required=True,
                    help="Input dir with chunks")
parser.add_argument("-u", "--chunk_filename_suffix", action="store", dest="chunk_filename_suffix", default=None,
                    help="Suffix of chunk filenames. "
                         "Chunks must be named as <prefix><separator><chunk number><suffix> . "
                         "Default: not set")
parser.add_argument("-p", "--chunk_filename_prefix", action="store", dest="chunk_filename_prefix", default=None,
                    help="Prefix of chunk filenames. "
                         "Chunks must be named as <prefix><separator><chunk number><suffix> . "
                         "Default: not set")
parser.add_argument("-s", "--separator", action="store", dest="separator", default="_",
                    help="Separator in chunk filename."
                         "Chunks must be named as <prefix><separator><chunk number><suffix> . "
                         "Default: '_'")
parser.add_argument("-n", "--total_number_of_chunks", action="store", dest="number_of_chunks", type=int,
                    required=True,
                    help="Total number of chunks")
parser.add_argument("-m", "--min_chunk_size", action="store", dest="min_chunk_size", type=int,
                    required=True,
                    help="Minimum size of chunk file.")

args = parser.parse_args()

AnnotationsRoutines.check_chunks(args.chunk_dir, args.number_of_chunks, args.min_chunk_size,
                                 separator=args.separator,
                                 chunk_filename_suffix=args.chunk_filename_suffix,
                                 chunk_filename_prefix=args.chunk_filename_prefix)
