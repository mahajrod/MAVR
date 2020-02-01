#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.Tools.WGA import LAST


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=sys.stdin,
                    help="Input MAF file")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output MAF file")

parser.add_argument("-e", "--max_e", action="store", dest="max_e", type=float,
                    help="Maximum e-value to keep alignment. Default: not set")
parser.add_argument("-g", "--max_eg2", action="store", dest="max_eg2", type=float,
                    help="Maximum eg2-value to keep alignment. Default: not set")
parser.add_argument("-s", "--min score", action="store", dest="min_score", type=float,
                    help="Minimum score to keep alignment. Default: not set")

args = parser.parse_args()

LAST.filter_lastall_maf(args.input, args.output, max_eg2=args.max_eg2, max_e=args.max_e, min_score=args.min_score)
