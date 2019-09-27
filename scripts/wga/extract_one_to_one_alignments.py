#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.WGA import LAST


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input MAF file")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files with")
parser.add_argument("-m", "--max_mismap_probability", action="store", dest="max_mismap_probability",
                    default=1.0, type=float,
                    help="Maximum mismatch probability. Default: 1.0")

args = parser.parse_args()


LAST.extract_one_to_one_alignments(args.input, args.output_prefix,
                                   max_mismap_probability=args.max_mismap_probability)
