#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Tools.Annotation import Barrnap

parser = argparse.ArgumentParser()

parser.add_argument("-r", "--report_file", action="store", dest="report_file", required=True,
                    help="File with barrnap report")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

args = parser.parse_args()


Barrnap.analyze_report(args.report_file, args.output_prefix)
