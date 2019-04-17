#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Tools.RepeatMasking import RepeatMasker

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input file with repeatmasking in gff format")
parser.add_argument("-p", "--output", action="store", dest="output", required=True,
                    help="Output file with repeat classes present in gff")
"""
parser.add_argument("-o", "--output_gff", action="store", dest="output_gff", required=True,
                    help="Output gff file")
parser.add_argument("-a", "--repeat_classes", action="store", dest="repeat_classes", required=True,
                    help="File to write annotated repeat classes")
parser.add_argument("-f", "--repeat_families", action="store", dest="repeat_families", required=True,
                    help="File to write annotated repeat families")
"""
args = parser.parse_args()

RepeatMasker.extract_annotated_repeat_types_from_gff(args.input, args.output)
