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
parser.add_argument("-l", "--label_list", action="store", dest="label_list", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of labels for species from MAF file")
parser.add_argument("-s", "--separator", action="store", dest="separator", default=".",
                    help="Separator between label and sequence id ")

args = parser.parse_args()

WGARoutines.label_maf(args.input, args.output, args.label_list, args.separator)
