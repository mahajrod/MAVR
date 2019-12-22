#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Parsers.CDHIT import CollectionCDHIT


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--clustering_file", action="store", dest="clustering_file", required=True,
                    help="Input file with CDHit clusters")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output fam file")
parser.add_argument("-f", "--format", action="store", dest="format", default="fam",
                    help="Format of output file. Allowed: fam, tab(default)")

args = parser.parse_args()

collection = CollectionCDHIT(input_file=args.clustering_file)
collection.write(args.output, format=args.format)
