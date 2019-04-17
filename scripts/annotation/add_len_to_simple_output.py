#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Annotation import Exonerate

parser = argparse.ArgumentParser()

parser.add_argument("-s", "--exonerate_simple_out", action="store", dest="exo_simple_out",
                    help="File with converted to 'simple' output")
parser.add_argument("-l", "--length_file", action="store", dest="len_file",
                    help="File with lengths of queues")
parser.add_argument("-o", "--output", action="store", dest="output",
                    help="Output file")

args = parser.parse_args()

Exonerate.add_len_to_simple_output(args.exo_simple_out, args.len_file, args.output)
