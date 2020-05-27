#!/usr/bin/env python
__author__ = 'mahajrod'

import argparse
from RouToolPa.Parsers.Bowtie2 import Bowtie2Table

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True, type=lambda s: s.split(","),
                    help="Comma-separated list of files with stats reporterd by Bowtie2")
parser.add_argument("-s", "--samples", action="store", dest="samples", required=True, type=lambda s: s.split(","),
                    help="Comma-separated list of samples")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

args = parser.parse_args()

"""

Usage:
 ~/Soft/MAVR/scripts/alignment/bowtie2/gather_stats_from_bowtie2_stat.py -i `ls ../assembly/*/*/*.fasta.stats | tr "\n" "," | sed s/,$//` -s `ls ../assembly  | tr "\n" "," | sed s/,$//` -o bowtie2.stats

"""
bowtie2_table = Bowtie2Table(args.input, samples=args.samples)
bowtie2_table.write("%s.tab" % args.output_prefix)
bowtie2_table.write_xlsx("%s.xlsx" % args.output_prefix)
bowtie2_table.draw(args.output_prefix, width=0.35, dpi=300, extensions=("png",))
