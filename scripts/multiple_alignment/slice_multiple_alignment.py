#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

from Bio import AlignIO


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input_file",
                    help="Input file with alignment")
parser.add_argument("-o", "--output_file", action="store", dest="output_file",
                    help="Output file with sliced alignment")
parser.add_argument("-s", "--start", action="store", dest="start", type=int,
                    help="Start. 1-based")
parser.add_argument("-e", "--end", action="store", dest="end", type=int,
                    help="End. 1-based")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of alignment")
args = parser.parse_args()

alignment = AlignIO.read(args.input_file, args.format)
AlignIO.write(alignment[:, args.start-1:args.end], args.output_file, args.format)
