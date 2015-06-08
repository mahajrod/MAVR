#!/usr/bin/env python
__author__ = 'mahajrod'
import argparse
import os
from Tools.Filter import TrimGalore
from Tools.Tuxedo import Tophat, Cufflinks

parser = argparse.ArgumentParser()

parser.add_argument("-l", "--left_reads", action="store", dest="left",
                    help="fastq file with left reads")

parser.add_argument("-r", "--right_reads", action="store", dest="right",
                    help="fastq file with right reads", default=None)

parser.add_argument("--left_trim", action="store", dest="left_trim",
                    help="5' trimm of left reads")

parser.add_argument("--right_trim", action="store", dest="right_trim",
                    help="5' trimm of right reads", default=None)

parser.add_argument("-a", "--adapter", action="store", dest="adapter",
                    help="adapter seqeunce for 3' filtering", default=None)

parser.add_argument("-q", "--quality", action="store", dest="quality",
                    help="base quality treshold, default=20", default=20)

parser.add_argument("-c", "--quality_score", action="store", dest="quality_score",
                    help="Type of quality score, allowed values: phred33, phred64", default="phred33")

parser.add_argument("-l", "--length", action="store", dest="length",
                    help="minimum length of reads after filtering", default=20)

parser.add_argument("-o", "--output", action="store", dest="output",
                    help="fasta file with reverse complement sequences")

args = parser.parse_args()

TrimGalore.filter(args.length, args.left, args.left_trim, args.right, args.right_trim, args.quality_score, args.adapter,
                  args.quality)

os.chdir("trimmed")
#Tophat.align()