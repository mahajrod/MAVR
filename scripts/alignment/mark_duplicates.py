#!/usr/bin/env python
__author__ = 'mahajrod'

import argparse
from RouToolPa.Tools.Picard import MarkDuplicates

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_bam", action="store", dest="input_bam", required=True,
                    help="Input BAM file")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Fasta file with query sequences")
parser.add_argument("-p", "--picard_dir", action="store", dest="picard_dir", default="",
                    help="Directory with PICARD jar")
parser.add_argument("-m", "--memory", action="store", dest="memory", default="100g",
                    help="Maximum memory to use. Default: 100g")
parser.add_argument("-t", "--tmp_dir", action="store", dest="tmp_dir",
                    help="Directory for temporary files. Default: system default")

args = parser.parse_args()

MarkDuplicates.path = args.picard_dir
MarkDuplicates.max_memory = args.memory
MarkDuplicates.tmp_dir = args.tmp_dir

MarkDuplicates.mkdup(args.input_bam, args.output_prefix)
