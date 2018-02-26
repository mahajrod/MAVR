#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse

from Tools.Alignment import LongRanger
from Routines import FileRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="File with input reference sequences")

parser.add_argument("-o", "--output", action="store", dest="output",
                    help="File to write prepared reference")

args = parser.parse_args()


LongRanger.prepare_reference(args.input, args.output, max_scaffold_length=527000000, max_scaffold_number=500,
                             polyN_len=500)
