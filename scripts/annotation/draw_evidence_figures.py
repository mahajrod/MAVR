#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

from Routines import DrawingRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input evidence file")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

args = parser.parse_args()

DrawingRoutines.draw_evidence_figures(args.input, args.output_prefix)
