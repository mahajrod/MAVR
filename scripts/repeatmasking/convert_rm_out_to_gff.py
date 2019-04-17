#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Tools.RepeatMasking import RepeatMasker

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input file with RepeatMasker output")
parser.add_argument("-p", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
"""
parser.add_argument("-o", "--output_gff", action="store", dest="output_gff", required=True,
                    help="Output gff file")
parser.add_argument("-a", "--repeat_classes", action="store", dest="repeat_classes", required=True,
                    help="File to write annotated repeat classes")
parser.add_argument("-f", "--repeat_families", action="store", dest="repeat_families", required=True,
                    help="File to write annotated repeat families")
"""
args = parser.parse_args()

output_gff = "%s.gff" % args.output_prefix
repeat_classes = "%s.repeat_classes" % args.output_prefix
repeat_families = "%s.repeat_families" % args.output_prefix

RepeatMasker.convert_rm_out_to_gff(args.input, output_gff, repeat_classes, repeat_families)
