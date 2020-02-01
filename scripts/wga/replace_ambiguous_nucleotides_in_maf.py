#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import WGARoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input MAF file")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output MAF file")


args = parser.parse_args()

WGARoutines.replace_ambigious_nucleotides_in_maf(args.input, args.output)
