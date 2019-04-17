#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Parsers.Primer3 import CollectionPrimer3

parser = argparse.ArgumentParser()

parser.add_argument("-o", "--primer3_output", action="store", dest="primer3_output", required=True,
                    help="File with Primer3 output")
parser.add_argument("-a", "--alignments_file", action="store", dest="alignments_file", required=True,
                    help="File to write alignments of primers to target sequence")

args = parser.parse_args()

primer3_collection = CollectionPrimer3(from_file=True, primer3_file=args.primer3_output)
primer3_collection.write_alignments(args.alignments_file, segment_length=120, left_primer_symbol=">",
                                    target_symbol="*", right_primer_symbol="<")

