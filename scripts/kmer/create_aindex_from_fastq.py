#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os
import argparse

from RouToolPa.Tools.AIndex import AIndex


parser = argparse.ArgumentParser()

parser.add_argument("-f", "--forward_file", action="store", dest="forward_file", required=True,
                    help="File with forward reads")
parser.add_argument("-r", "--reverse_file", action="store", dest="reverse_file", required=True,
                    help="File with reverse reads")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Output prefix")
parser.add_argument("-m", "--kmer_length", action="store", dest="kmer_length", type=int, default=23,
                    help="Length of kmers. Default - 23")
parser.add_argument("-s", "--hash_size", action="store", dest="hash_size", type=str, default="1000000",
                    help="Size of hash. Estimation of hash size: for short reads S=(G + k*n)/0.8, "
                    "G - genome size, k - kmer length, n - number of reads, for assembled sequences "
                    "S=Sum(L). Default - 1000000")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads. Default - 1")

args = parser.parse_args()

AIndex.threads = args.threads

AIndex.create_index_from_fastq(args.forward_file, args.reverse_file,
                               args.kmer_length, args.output_prefix,
                               lower_count=2, upper_count=None,
                               filetype="fastq",
                               create_aindex=True,
                               hash_size=args.hash_size)
