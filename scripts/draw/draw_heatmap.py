#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import MatplotlibRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input_file",
                    help="Input file with data")
parser.add_argument("-c", "--x_column", action="store", dest="x_col", default=0, type=int,
                    help="Column with x values. Default - 0")
parser.add_argument("-d", "--y_column", action="store", dest="y_col", default=1, type=int,
                    help="Column with y values. Default - 1")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    help="Prefix of output files")
#parser.add_argument("-s", "--separator", action="store", dest="separator", default="\n",
#                    help="Separator between values in input file. Default - '\\n', i.e. one value per line")
parser.add_argument("-b", "--number_of_bins", action="store", dest="number_of_bins",
                    type=lambda s: map(int, s.split(",") if "," in s else int(s)), default=20,
                    help="Number of bins in histogram. Can be either integer or comma-separated list of two integers."
                         "Incompatible with -w/--width_of_bins option. Default - 20")
parser.add_argument("-w", "--width_of_bins", action="store", dest="width_of_bins",
                    type=lambda s: map(int, s.split(",") if "," in s else int(s)),
                    help="Width of bins in histogram. Can be either integer or comma-separated list of two integers."
                         "Incompatible with -a/--array_of_bins option. Not set by default")
parser.add_argument("-a", "--array_of_bins", action="store", dest="array_of_bins",
                    type=lambda s: map(lambda b: b.split(","), s.split(";")) if ";" in s else s.split(","),
                    help="Array of bins in histogram. Can be either comma-separate list of integers or "
                         "two comma-separated lists separated by dot with comma"
                         "Incompatible with -w/--width_of_bins option. Not set by default")
parser.add_argument("-m", "--min_counts_to_show", action="store", dest="min_counts_to_show", type=int, default=1,
                    help="Minimum value to show. Default - 1")
parser.add_argument("-s", "--remove_colorbar", action="store_true", dest="remove_colorbar", default=False,
                    help="Remove colorbar. Default - False")
parser.add_argument("-n", "--min_x", action="store", dest="min_x", type=float,
                    help="Minimum x value to show. Default - min(x)")
parser.add_argument("-x", "--max_x", action="store", dest="max_x", type=float,
                    help="Maximum x value to show. Default - max(x)")

parser.add_argument("-z", "--min_y", action="store", dest="min_y", type=float,
                    help="Minimum y value to show. Default - min(y)")
parser.add_argument("-q", "--max_y", action="store", dest="max_y", type=float,
                    help="Maximum y value to show. Default - max(y)")
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

MatplotlibRoutines.draw_heatmap_from_file(args.input_file, args.output_prefix, x_column=args.x_col, y_column=args.y_col,
                                          xlabel=args.xlabel, ylabel=args.ylabel, title=args.title,
                                          figsize=(8, 8), minimum_counts_to_show=args.min_counts_to_show,
                                          extensions=args.extensions, show_colorbar=not args.remove_colorbar,
                                          bin_number=args.number_of_bins, bin_width=args.width_of_bins,
                                          bin_array=args.array_of_bins,
                                          min_x_value=args.min_x, max_x_value=args.max_x,
                                          min_y_value=args.min_y, max_y_value=args.max_y,
                                          add_max_value=True)

