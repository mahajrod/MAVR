#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse

from Tools.Annotation import Exonerate
from Routines.File import make_list_of_path_to_files

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    type=lambda s: make_list_of_path_to_files(s.split(",")),
                    help="Input comma-separated list of files/directories with exonerate output")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    help="Prefix of output files")

args = parser.parse_args()

Exonerate.split_output(args.input, args.output_prefix)
