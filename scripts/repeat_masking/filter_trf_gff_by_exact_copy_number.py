#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import shutil
import argparse

from Tools.RepeatMasking import TRF

from Routines.File import split_filename


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input gff with TRF repeats")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output gff with filtered repeats")
parser.add_argument("-x", "--filtered_out", action="store", dest="filtered_out", required=True,
                    help="Gff with filtered out repeats")
parser.add_argument("-b", "--min_copy_number", action="store", dest="min_copy_number", type=int, required=True,
                    help="Minimum number of consequent exact copies to extract")

args = parser.parse_args()


TRF.filter_trf_gff_by_exact_copy_number(args.input, args.output, args.filtered_out, args.min_copy_number)

