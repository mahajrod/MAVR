#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
#import os
import argparse
"""
import numpy as np

import matplotlib
matplotlib.use('Agg')
os.environ['MPLCONFIGDIR'] = '/tmp/'
import matplotlib.pyplot as plt

from Bio import SeqIO

from RouToolPa.Routines import MatplotlibRoutines
from RouToolPa.Routines.Sequence import rev_com_generator
from RouToolPa.Routines.File import make_list_of_path_to_files
"""
from RouToolPa.Tools.Kmers import Jellyfish


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with data for histogram")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Output prefix")
parser.add_argument("-e", "--output_formats", action="store", dest="output_formats", type=lambda s: s.split(","),
                    default=["svg", "png"],
                    help="Comma-separated list of formats (supported by matlotlib) "
                         "of output figure.Default: svg,eps,pdf,png,jpg")
parser.add_argument("-l", "--logbase", action="store", dest="logbase", type=int, default=10,
                    help="Base of logarithm. Default - 10")
parser.add_argument("-w", "--low_limit", action="store", dest="low_limit", type=int, default=5,
                    help="Low limit of histogram without logscale")
parser.add_argument("-g", "-high_limit", action="store", dest="high_limit", type=int, default=100,
                    help="High limit of histogram without logscale")
parser.add_argument("-m", "--kmer_length", action="store", dest="kmer_length", type=int, default=23,
                    help="Length of kmers. Default - 23")

#parser.add_argument("-d", "--draw_peaks_and_gaps", action="store_true", dest="draw_peaks_and_gaps",
#                    help="Draw peaks and gaps")

args = parser.parse_args()

Jellyfish.draw_kmer_distribution(args.input, args.kmer_length, args.output_prefix, output_formats=args.output_formats,
                                 logbase=args.logbase, non_log_low_limit=args.low_limit,
                                 non_log_high_limit=args.high_limit) #, draw_peaks_and_gaps=args.draw_peaks_and_gaps)
"""
bins, counts = np.loadtxt(args.input, unpack=True)

figure = plt.figure(1, figsize=(8, 8), dpi=300)
subplot = plt.subplot(1, 1, 1)
plt.suptitle("Distribution of %i-mers" % args.kmer_length, fontweight='bold')
plt.plot(bins, counts)
plt.xlim(xmin=1, xmax=10000000)
plt.xlabel("Multiplicity")
plt.ylabel("Number of distinct %s-mers" % args.kmer_length)
subplot.set_yscale('log', basey=args.logbase)
subplot.set_xscale('log', basex=args.logbase)

for extension in args.output_formats:
    plt.savefig("%s.logscale.%s" % (args.output_prefix, extension))

plt.close()

selected_counts = counts[args.low_limit-1:args.high_limit]
selected_bins = bins[args.low_limit-1:args.high_limit]

figure = plt.figure(2, figsize=(8, 8), dpi=300)
subplot = plt.subplot(1, 1, 1)
plt.suptitle("Distribution of %s-mers" % args.kmer_length, fontweight='bold')
plt.plot(selected_bins, selected_counts)

plt.xlabel("Multiplicity")
plt.ylabel("Number of distinct %s-mers" % args.kmer_length)
plt.xlim(xmin=args.low_limit, xmax=args.high_limit)

for extension in args.output_formats:
    plt.savefig("%s.no_logscale.%s" % (args.output_prefix, extension))

plt.close()

figure = plt.figure(3, figsize=(6, 12), dpi=400)
subplot_list = []
for i, b, c in zip([1, 2], [bins, selected_bins], [counts, selected_counts]):
    subplot_list.append(plt.subplot(2, 1, i))
    plt.suptitle("Distribution of %s-mers" % args.kmer_length, fontweight='bold', fontsize=13)
    plt.plot(b, c)

    plt.ylabel("Number of distinct %s-mers" % args.kmer_length, fontsize=13)
    if i == 1:
        subplot_list[0].set_yscale('log', basey=args.logbase)
        subplot_list[0].set_xscale('log', basex=args.logbase)
        plt.xlim(xmin=1, xmax=10000000)
    elif i == 2:
        plt.xlim(xmin=args.low_limit, xmax=args.high_limit)
        plt.xlabel("Multiplicity", fontsize=15)

MatplotlibRoutines.zoom_effect(subplot_list[0], subplot_list[1], args.low_limit, args.high_limit)
plt.subplots_adjust(hspace=0.12, wspace=0.05, top=0.95, bottom=0.05, left=0.14, right=0.95)

for extension in args.output_formats:
    plt.savefig("%s.%s" % (args.output_prefix, extension))
"""
