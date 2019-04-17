#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Expression import Gffread


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input gff with annotation")
parser.add_argument("-g", "--genome", action="store", dest="genome", required=True,
                    help="Fasta with genome")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file to write cds")

args = parser.parse_args()

Gffread.extract_cds(args.input, args.genome, args.output)
