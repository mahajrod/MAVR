#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import GORoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="File with annotations made by EggNOG")
parser.add_argument("-g", "--go_file", action="store", dest="go_file", required=True,
                    help="File with GO terms")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

args = parser.parse_args()


GORoutines.extract_entries_by_GO_from_eggnogmapper_output(args.input, args.go_file, args.output_prefix,
                                                          comments_prefix="#", separator="\t")

