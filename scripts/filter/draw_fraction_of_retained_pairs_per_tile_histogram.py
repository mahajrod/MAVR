#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Parsers.FaCut import FaCutReport


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--facut_report", action="store", dest="facut_report", required=True,
                    help="File with FaCut report")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    default="./", help="Prefix of output files")

args = parser.parse_args()

facut_report = FaCutReport(args.facut_report)

facut_report.draw_fraction_of_retained_pairs_per_tile_histogram(args.output_prefix)
