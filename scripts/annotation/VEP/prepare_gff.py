#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Annotation import VEP

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_gff", action="store", dest="input", required=True, type=lambda s: s.split(","),
                    help="List of input .gff file")
parser.add_argument("-p", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

args = parser.parse_args()

VEP.prepare_gff(args.input, args.output_prefix)
