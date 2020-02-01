#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.Routines import WGARoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input MAF file")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output MAF file")


args = parser.parse_args()

WGARoutines.replace_ambigious_nucleotides_in_maf(args.input, args.output)
