#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Parsers.GFF import CollectionGFF

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input bed file")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="File to write GFF output")

args = parser.parse_args()

bed_data = CollectionGFF(in_file=args.input, parsing_mode="only_coordinates", format="bed")

bed_data.write(args.output, output_format="gff")
