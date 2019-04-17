#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import SequenceRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with sequences")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file with corrected sequences")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of input and output files. Allowed formats genbank, fasta(default)")
parser.add_argument("-b", "--black_list", action="store", dest="black_list",
                    help="File with record ids from black list")
parser.add_argument("-w", "--white_list", action="store", dest="white_list",
                    help="File with record ids from white list")
parser.add_argument("-m", "--masking", action="store", dest="masking",
                    help="0-based BED file with regions to mask")
parser.add_argument("-t", "--trimming", action="store", dest="trimming",
                    help="0-based BED file with regions to trim")
args = parser.parse_args()

SequenceRoutines.correct_sequences_from_file(args.input, args.output,
                                             black_list_file=args.black_list,
                                             white_list_file=args.white_list,
                                             regions_to_trim_file=args.trimming,
                                             regions_to_mask_file=args.masking,
                                             parsing_mode="parse",
                                             format=args.format)
