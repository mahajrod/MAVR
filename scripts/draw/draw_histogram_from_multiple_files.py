#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import DrawingRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input",
                    type=lambda s: DrawingRoutines.make_list_of_path_to_files(s.split(",")),
                    help="Comma-separated list of files and directories")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    help="Prefix of output files")
parser.add_argument("-s", "--separator", action="store", dest="separator", default="\n",
                    help="Separator between values in input file. Default - '\\n', i.e. one value per line")

parser.add_argument("-b", "--number_of_bins", action="store", dest="number_of_bins", type=int, default=30,
                    help="Number of bins in histogram. Incompatible with -w/--width_of_bins option. Default - 30")
"""
parser.add_argument("-w", "--width_of_bins", action="store", dest="width_of_bins", type=float,
                    help="Width of bins in histogram. Incompatible with -b/--number_of_bins option. Not set by default")
"""
parser.add_argument("-n", "--min_value", action="store", dest="min", type=float, default=0,
                    help="Minimum value to show. Default - 1")
parser.add_argument("-x", "--max_value", action="store", dest="max", type=float,
                    help="Maximum value to show. Default - length of longest sequence")

"""
parser.add_argument("-g", "--logbase", action="store", dest="logbase", type=int, default=10,
                    help="Logbase to use for log-scaled histograms")
"""

parser.add_argument("-e", "--extensions", action="store", dest="extensions", type=lambda x: x.split(","),
                    default=["png"],
                    help="Comma-separated list of extensions for histogram files. Default: png only")
parser.add_argument("-l", "--xlabel", action="store", dest="xlabel",
                    help="X label")
parser.add_argument("-y", "--ylabel", action="store", dest="ylabel",
                    help="Y label")
parser.add_argument("-t", "--title", action="store", dest="title",
                    help="Title of histogram")

args = parser.parse_args()

DrawingRoutines.draw_histogram_from_multiple_files(args.input, args.output_prefix,  filelabels=None,
                                                   nbins=args.number_of_bins,
                                                   figsize=(5, 5), title=args.title, xlabel=args.xlabel,
                                                   ylabel=args.ylabel,
                                                   extensions=args.extensions,
                                                   separator=args.separator,
                                                   xmin=args.min, xmax=args.max,
                                                   histtype="stepfilled")
