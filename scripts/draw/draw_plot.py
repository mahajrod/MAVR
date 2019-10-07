#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

import numpy as np
import scipy.stats as stats

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input_file",
                    help="Input file with data")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    help="Prefix of output files")
parser.add_argument("-s", "--separator", action="store", dest="separator", default="\t",
                    help="Separator between values in input file. Default - '\\t' ")

parser.add_argument("-a", "--x_column_index", action="store", dest="x_column_index", type=int,
                    help="Index of column with x values. 0-based")
parser.add_argument("-b", "--y_column_index", action="store", dest="y_column_index", type=int,
                    help="Index of column with y values. 0-based")

parser.add_argument("-n", "--min_x", action="store", dest="min_x", type=float,
                    help="Minimum x value to show. Default - not set")
parser.add_argument("-x", "--max_x", action="store", dest="max_x", type=float,
                    help="Maximum x value to show. Default - not set")
parser.add_argument("-q", "--min_y", action="store", dest="min_y", type=float,
                    help="Minimum y value to show. Default - not set")
parser.add_argument("-r", "--max_y", action="store", dest="max_y", type=float,
                    help="Maximum y value to show. Default - not set")
parser.add_argument("-e", "--extensions", action="store", dest="extensions", type=lambda x: x.split(","),
                    default=["png", "svg"],
                    help="Comma-separated list of extensions for histogram files")
parser.add_argument("-l", "--xlabel", action="store", dest="xlabel",
                    help="X label")
parser.add_argument("-y", "--ylabel", action="store", dest="ylabel",
                    help="Y label")
parser.add_argument("-t", "--title", action="store", dest="title",
                    help="Title of histogram")
parser.add_argument("--width", action="store", dest="width", default=6, type=int,
                    help="Figure width. Default: 6")
parser.add_argument("--height", action="store", dest="height", default=6, type=int,
                    help="Figure height. Default: 6")
parser.add_argument("-m", "--markersize", action="store", dest="markersize", default=2, type=int,
                    help="Size of marker. Default: 2")
parser.add_argument("--ylog", action="store", dest="ylogbase", default=10, type=int,
                    help="Log base for figure with logarithmic scale on y axis. Default: 10")
parser.add_argument("--type", action="store", dest="type", default="plot",
                    help="Type of figure. Allowed: plot(default), scatter")
parser.add_argument("-g", "--grid", action="store_true", dest="grid",
                    help="Show grid. Default: False")
args = parser.parse_args()

data = np.loadtxt(args.input_file, comments="#", usecols=(args.x_column_index, args.y_column_index))

plt.figure(1, figsize=(args.width, args.height), dpi=300)
plt.subplot(1, 1, 1)
if args.type == "plot":
    plt.plot(data[:, 0], data[:, 1], markersize=args.markersize)
elif args.type == "scatter":
    plt.scatter(data[:, 0], data[:, 1], s=args.markersize)
plt.xlim(xmin=args.min_x, xmax=args.max_x)
plt.ylim(ymin=args.min_y, ymax=args.max_y)
if args.xlabel:
    plt.xlabel(args.xlabel)
if args.ylabel:
    plt.ylabel(args.ylabel)
if args.title:
    plt.title(args.title)
if args.grid:
    plt.grid()
print("Kendal's tau")
print(stats.kendalltau(data[:, 0], data[:, 1]))

print("Pearson's r")
print(stats.pearsonr(data[:, 0], data[:, 1]))
for ext in args.extensions:
    plt.savefig("%s.%s.%s" % (args.output_prefix, args.type, ext))

plt.yscale("log")

for ext in args.extensions:
    plt.savefig("%s.%s.ylog%i.%s" % (args.output_prefix, args.type, args.ylogbase, ext))
