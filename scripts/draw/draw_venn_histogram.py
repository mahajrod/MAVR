#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()

from RouToolPa.Routines import MatplotlibRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-a", "--id_file_a", action="store", dest="id_file_a", required=True,
                    help="Id File A ")
parser.add_argument("-b", "--id_file_b", action="store", dest="id_file_b", required=True,
                    help="Id File B ")
parser.add_argument("-c", "--id_file_c", action="store", dest="id_file_c",
                    help="Id File C ")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-l", "--set_labels", action="store", dest="set_labels", type=lambda s: s.split(","),
                    help="Comma-separated list of set labels")
parser.add_argument("-r", "--set_colors", action="store", dest="set_colors", type=lambda s: s.split(","),
                    help="Comma-separated list of set colors")

"""
parser.add_argument("-s", "--separator", action="store", dest="separator", default="\t",
                    help="Separator between values in input file. Default - '\\t' ")

parser.add_argument("-a", "--x_column_index", action="store", dest="x_column_index", type=int,
                    help="Index of column with x values. 0-based")
parser.add_argument("-b", "--y_column_index", action="store", dest="y_column_index", type=int,
                    help="Index of column with y values. 0-based")

parser.add_argument("-n", "--min_value", action="store", dest="min_length", type=float, default=0,
                    help="Minimum value to show. Default - 1")
parser.add_argument("-x", "--max_value", action="store", dest="max_length", type=float,
                    help="Maximum value to show. Default - length of longest sequence")
"""
parser.add_argument("-e", "--extensions", action="store", dest="extensions", type=lambda x: x.split(","),
                    default=["png", "svg"],
                    help="Comma-separated list of extensions for histogram files")
"""
parser.add_argument("-l", "--xlabel", action="store", dest="xlabel",
                    help="X label")
parser.add_argument("-y", "--ylabel", action="store", dest="ylabel",
                    help="Y label")
"""
parser.add_argument("-t", "--title", action="store", dest="title",
                    help="Title of histogram")

args = parser.parse_args()

MatplotlibRoutines.venn_diagram_from_sets_from_files(args.id_file_a, args.id_file_b, set3_file=args.id_file_c,
                                                     set_labels=args.set_labels, set_colors=args.set_colors,
                                                     output_prefix=args.output_prefix, extensions=args.extensions,
                                                     title=args.title)
