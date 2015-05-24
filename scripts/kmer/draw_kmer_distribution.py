#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os
import argparse

import matplotlib.pyplot as plt

from Tools.Kmers import Jellyfish

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input file with kmer database.")
parser.add_argument("-o", "--output_prefix", action="store", dest="out_prefix",
                    help="Prefix_of_output_files")
parser.add_argument("-e", "--output_formats", action="store", dest="output_formats", type=lambda s: s.split(","),
                    default=["svg", "eps", "pdf", "png", "jpg"],
                    help="Comma-separated list of formats (supported by matlotlib) "
                         "of output figure.Default: svg,eps,pdf,png,jpg")
parser.add_argument("-w", "--low_limit", action="store", dest="low_limit", type=int, default=1,
                    help="Low limit of histogram")
parser.add_argument("-g", "-high_limit", action="store", dest="high_limit", type=int, default=10000,
                    help="High limit of histogram")
parser.add_argument("-b", "--bin_width", action="store", dest="bin_width", type=int, default=1,
                    help="Bin width of histogram")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads")
parser.add_argument("-d", "--_without_logscale", action="store_true", dest="without_logscale",
                    help="Dont logsscale axes")
parser.add_argument("-a", "--logbase", action="store", dest="logbase", type=int, default=10,
                    help="Base of logarithm")
parser.add_argument("-k", "--kmer_size", action="store", dest="kmer_size", type=int,
                    help="Size of kmers in base. Used in suptitle of figure.")
args = parser.parse_args()

file_prefix = ".".join(os.path.basename(args.input).split(".")[:-1])
histo_file = "%s.histo" % file_prefix

if args.out_prefix is None:
    args.out_prefix = file_prefix + "_histogram"

Jellyfish.threads = args.threads
Jellyfish.histo(args.input, histo_file, bin_width=args.bin_width,
                lower_count=args.low_limit, upper_count=args.high_limit)

counts = []
bins = []

with open(histo_file, "r") as histo_fd:
    for line in histo_fd:
        entry = line.strip().split()
        counts.append(entry[1])
        bins.append(entry[0])


figure = plt.figure(1, figsize=(8, 8), dpi=300)
subplot = plt.subplot(1, 1, 1)
plt.suptitle("Distribution of %s-mers" % (str(args.kmer_size) if args.kmer_size else "k"),
             fontweight='bold')
plt.plot(bins, counts)

plt.xlabel("Multiplicity")
plt.ylabel("Number of distinct kmers")
if not args.without_logscale:
    subplot.set_yscale('log', basey=args.logbase)
    subplot.set_xscale('log', basex=args.logbase)

for extension in args.output_formats:
    plt.savefig("%s.%s" % (args.out_prefix, extension))



