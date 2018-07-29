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

parser.add_argument("-n", "--min_value", action="store", dest="min_length", type=float, default=0,
                    help="Minimum value to show. Default - 1")
parser.add_argument("-x", "--max_value", action="store", dest="max_length", type=float,
                    help="Maximum value to show. Default - length of longest sequence")
parser.add_argument("-e", "--extensions", action="store", dest="extensions", type=lambda x: x.split(","),
                    default=["png", "svg"],
                    help="Comma-separated list of extensions for histogram files")
parser.add_argument("-l", "--xlabel", action="store", dest="xlabel",
                    help="X label")
parser.add_argument("-y", "--ylabel", action="store", dest="ylabel",
                    help="Y label")
parser.add_argument("-t", "--title", action="store", dest="title",
                    help="Title of histogram")

args = parser.parse_args()

data = np.loadtxt(args.input_file, comments="#")
print

plt.figure(1, figsize=(6, 6))
plt.subplot(1, 1, 1)

plt.plot(data[:, args.x_column_index], data[:, args.y_column_index], "b.")
plt.xlim(xmin=args.min_length, xmax=args.max_length)
if args.xlabel:
    plt.xlabel(args.xlabel)
if args.ylabel:
    plt.ylabel(args.ylabel)
if args.title:
    plt.title(args.title)

print("Kendal's tau")
print(stats.kendalltau(data[:, args.x_column_index], data[:, args.y_column_index]))

print("Pearson's r")
print(stats.pearsonr(data[:, args.x_column_index], data[:, args.y_column_index]))
for ext in args.extensions:
    plt.savefig("%s.%s" % (args.output_prefix, ext))

