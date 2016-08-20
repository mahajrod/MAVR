#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with sequences")
parser.add_argument("-s", "--separator", action="store", dest="separator", default="@",
                    help="Separator to use. default - '@'")
parser.add_argument("-o", "--output", action="store", dest="output", default="stdout",
                    help="File to write sequences with truncated ids. Default - stdout")

args = parser.parse_args()

out_file = "" if args.output == "stdout" else " > %s" % args.output

sed_string = "sed 's/^>\\(.*\\)%s.*/>\\1/' %s %s" % (args.separator, args.input, out_file)

os.system(sed_string)
