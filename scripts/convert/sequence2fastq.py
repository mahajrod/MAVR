#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from BConverters import SequenceConverters

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input sequence file")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output fastq file")
parser.add_argument("-d", "--default_quality", action="store", dest="default_quality", required=True, type=int,
                    help="Quality for bases in output fastq file")
parser.add_argument("-m", "--mode", action="store", dest="mode", default="index_db",
                    help="Parsing mode for sequence file. Allowed: parse, index, index_db(default)")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of sequence file. Default: fasta")
parser.add_argument("--index_file", action="store", dest="index_file",
                    help="Index file for parsing in index_db mode. Default: construct new")
parser.add_argument("-q", "--quality_score", action="store", dest="quality_score", default=33, type=int,
                    help="Quality score(offset). Allowed: 33(default), 64")

args = parser.parse_args()

SequenceConverters.sequence2fastq(args.input, args.output, args.mode, args.default_quality, format=args.format,
                                  index_file=args.index_file, score_type=args.quality_score)
