#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.BLAST import BLASTPlus

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input with blast hits in table format")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file")

args = parser.parse_args()

BLASTPlus.extract_hits_from_tbl_output(args.input, args.output)
