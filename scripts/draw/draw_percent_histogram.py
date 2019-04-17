#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import MatplotlibRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input_file",
                    help="Input file with data")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    help="Prefix of output files")
parser.add_argument("-s", "--separator", action="store", dest="separator",
                    help="Separator between values in input file. Default any whitespace")
parser.add_argument("-d", "--data_type", action="store", dest="data_type",
                    help="Data type. Default float")
parser.add_argument("-c", "--comments_prefix", action="store", dest="comments_prefix", default="#",
                    help="Prefix of lines with comments. Default - '#'")
parser.add_argument("-b", "--number_of_bins", action="store", dest="number_of_bins", type=int, default=20,
                    help="Number of bins in histogram. Default - 20")
parser.add_argument("-l", "--columns_list", action="store", dest="columns_list", type=lambda x: map(int, x.split(",")),
                    help="Comma-separated list of columns(0-based) with data. Default - all")
parser.add_argument("-e", "--extensions", action="store", dest="extensions", type=lambda x: x.split(","),
                    default=["png", "svg"],
                    help="Comma-separated list of extensions for histogram files")
parser.add_argument("-x", "--xlabel", action="store", dest="xlabel", default="%%",
                    help="X label. Default - '%%' ")
parser.add_argument("-y", "--ylabel", action="store", dest="ylabel", default="Number",
                    help="Y label. Default - 'Number'")
parser.add_argument("-t", "--title", action="store", dest="title",
                    help="Title of histogram")
parser.add_argument("-a", "--legend_location", action="store", dest="legend_location", default="upper center",
                    help="Location of legend on histogram. Default - 'upper center'")
parser.add_argument("-m", "--input_mode", action="store", dest="input_mode", default="percent",
                    help="Type of input data. Allowed: fraction, percent. Default - percent")
args = parser.parse_args()

MatplotlibRoutines.percent_histogram_from_file(args.input_file, args.output_prefix, data_type=args.data_type,
                                               column_list=args.columns_list, separator=args.separator,
                                               comments=args.comments_prefix, n_bins=args.number_of_bins,
                                               title=args.title, xlabel=args.xlabel, ylabel=args.ylabel,
                                               extensions=args.extensions, legend_location=args.legend_location,
                                               stats_as_legend=True, input_mode=args.input_mode)
