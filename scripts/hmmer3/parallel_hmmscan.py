#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse

from Tools.HMMER import HMMER3

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_hmm", action="store", dest="input",
                    help="Input hmm3 file")
parser.add_argument("-s", "--input_seq", action="store", dest="input_seq",
                    help="Input file with sequences")
parser.add_argument("-o", "--output_file", action="store", dest="output",
                    help="Output file")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads")
args = parser.parse_args()


HMMER3.threads = 1
HMMER3.parallel_hmmscan(args.input, args.input_seq, args.output, num_of_seqs_per_scan=None, split_dir="splited_fasta",
                        splited_output_dir="splited_output_dir", threads=args.threads)