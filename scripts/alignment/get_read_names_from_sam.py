#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Tools.Samtools import SamtoolsV1



parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input",
                    help="Input sam file")
parser.add_argument("-o", "--output", action="store", dest="output",
                    help="Output file with read names")

args = parser.parse_args()

SamtoolsV1.get_read_names(args.input, args.output)
