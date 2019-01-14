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

parser.add_argument("-a", "--column_index_list", action="store", dest="column_index_list",
                    type=lambda x: map(int, x.split(",")),
                    help="Comma-separated list of indexes for columns with values. 0-based")

parser.add_argument("-n", "--min_value", action="store", dest="min_value", type=float,
                    help="Minimum value to show. Default: not applied")
parser.add_argument("-x", "--max_value", action="store", dest="max_value", type=float,
                    help="Maximum value to show. Default: not applied")
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

data = np.loadtxt(args.input_file, comments="#", usecols=args.column_index_list)

if args.min_value and args.max_value:
    data = data[(data <= args.min_value) & (data >= args.max_value)]
elif args.max_value:
    data = data[data <= args.max_value]
elif args.min_value:
    data = data[data >= args.min_value]

plt.figure(1, figsize=(6, 6))
plt.subplot(1, 1, 1)
print data
print len(data)
plt.bar(np.arange(1, len(data) + 1, 1), data)
plt.xlim(xmin=0, xmax=len(data))

if args.xlabel:
    plt.xlabel(args.xlabel)
if args.ylabel:
    plt.ylabel(args.ylabel)
if args.title:
    plt.title(args.title)

for ext in args.extensions:
    plt.savefig("%s.%s" % (args.output_prefix, ext))
