#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

from Tools.RepeatMasking import RepeatMasker

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input file with RepeatMasker output")
parser.add_argument("-o", "--output_gff", action="store", dest="output_gff", required=True,
                    help="Output gff file")
parser.add_argument("-a", "--repeat_classes", action="store", dest="repeat_classes", required=True,
                    help="File to write annotated repeat classes")
parser.add_argument("-f", "--repeat_families", action="store", dest="repeat_families", required=True,
                    help="File to write annotated repeat families")

args = parser.parse_args()

RepeatMasker.convert_rm_out_to_gff(args.input, args.output_gff, args.repeat_classes, args.repeat_families)
