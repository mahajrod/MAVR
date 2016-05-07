#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse

from Routines import TreeFamRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with proteins")
parser.add_argument("-s", "--species_name", action="store", dest="species_name", required=True,
                    help="Species name to use as label")
parser.add_argument("-o", "--output", action="store", dest="output", default="stdout",
                    help="File to labeled proteins. Default - stdout")

args = parser.parse_args()

out_file = "" if args.output == "stdout" else " > %s" % args.output

sed_string = "sed 's/^>/>%s./' %s %s" % (args.species_name, args.input, out_file)
