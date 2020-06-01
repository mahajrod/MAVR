#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with sequences")
parser.add_argument("-l", "--label", action="store", dest="label", required=True,
                    help="Label to use")
parser.add_argument("-s", "--separator", action="store", dest="separator", default="@",
                    help="Separator to use. default - '@'")
parser.add_argument("-o", "--output", action="store", dest="output", default="stdout",
                    help="File to write unlabeled sequences. Default - stdout")

args = parser.parse_args()

out_file = "" if args.output == "stdout" else " > %s" % args.output

if args.input[-3:] == ".gz":
    sed_string = "zcat %s | sed 's/^>.*%s/>/'  %s" % (args.input, args.separator, out_file)
else:
    sed_string = "sed 's/^>.*%s/>/' %s %s" % (args.separator, args.input, out_file)

os.system(sed_string)
