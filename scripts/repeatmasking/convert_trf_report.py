#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Tools.RepeatMasking import TRF

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with TRF report")
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

TRF.convert_trf_report(args.input, args.output_prefix)
