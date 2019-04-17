#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from collections import OrderedDict
from RouToolPa.Parsers.KrATER import KrATERReportCollection

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of files with KrATER reports")
parser.add_argument("-s", "--samples", action="store", dest="samples",
                    type=lambda s: s.split(","),
                    help="Comma-separated list of sample ids. Must have same length as list of files."
                         "If not set filenames will be treated as sample ids")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="File to write general statistics")

args = parser.parse_args()

filename_dict = OrderedDict()
if args.samples:
    for i in range(0, len(args.input)):
        filename_dict[args.samples[i]] = args.input[i]
else:
    for filename in args.input:
        filename_dict[filename] = filename

report_collection = KrATERReportCollection(file_dict=filename_dict)

report_collection.write_general_stats(args.output)
