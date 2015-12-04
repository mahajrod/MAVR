#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse

from Tools.Annotation import Exonerate



parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with exonerate output")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    help="Prefix of output files")

args = parser.parse_args()

Exonerate.split_output(args.input, args.output_prefix)
