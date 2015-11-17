#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

from Tools.RepeatMasking import RepeatMasker

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input file with RepeatMasker output")

parser.add_argument("-o", "--output_gff", action="store", dest="output_gff", required=True,
                    help="Output gff file")

args = parser.parse_args()

RepeatMasker.convert_rm_out_to_gff(args.input, args.output_gff)


