#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", default="stdin",
                    help="Input file. Default: stdin")
parser.add_argument("-o", "--output_file", action="store", dest="output", default="stdout",
                    help="Output file. Default: stdout")

args = parser.parse_args()

in_fd = sys.stdin if args.input == "stdin" else open(args.input, "r")
out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

if args.output != "output":
    out_fd.close()
if args.input != "stdin":
    in_fd.close()