#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse

from RouToolPa.Routines import FileRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file")
parser.add_argument("-c", "--column_index", action="store", dest="column_index", required=True, type=int,
                    help="Column to label")
parser.add_argument("-l", "--label", action="store", dest="label", required=True,
                    help="Label to use")
parser.add_argument("-s", "--column_separator", action="store", dest="column_separator", default="\t",
                    help="Separator used in input file. Default: TAB")
parser.add_argument("-a", "--label_separator", action="store", dest="label_separator", default="@",
                    help="Label separator. Default: @")
parser.add_argument("-p", "--label_position", action="store", dest="label_position", default="first",
                    help="Label position. Allowed: first(default), last")

args = parser.parse_args()

FileRoutines.label_column_in_file(args.input, args.label, args.column_index, args.output,
                                  column_separator=args.column_separator,
                                  label_position=args.label_position,
                                  label_separator=args.label_separator,
                                  comments_prefix="#")
