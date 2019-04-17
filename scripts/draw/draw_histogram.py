#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import MatplotlibRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input_file",
                    help="Input file with data")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    help="Prefix of output files")
parser.add_argument("-s", "--separator", action="store", dest="separator", default="\n",
                    help="Separator between values in input file. Default - '\\n', i.e. one value per line")

parser.add_argument("-b", "--number_of_bins", action="store", dest="number_of_bins", type=int,
                    help="Number of bins in histogram. Incompatible with -w/--width_of_bins option. Default - 30")
parser.add_argument("-w", "--width_of_bins", action="store", dest="width_of_bins", type=float,
                    help="Width of bins in histogram. Incompatible with -b/--number_of_bins option. Not set by default")
parser.add_argument("-n", "--min_value", action="store", dest="min_length", type=float, default=0,
                    help="Minimum value to show. Default - 1")
parser.add_argument("-x", "--max_value", action="store", dest="max_length", type=float,
                    help="Maximum value to show. Default - length of longest sequence")
parser.add_argument("-g", "--logbase", action="store", dest="logbase", type=int, default=10,
                    help="Logbase to use for log-scaled histograms")
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

MatplotlibRoutines.draw_histogram_from_file(args.input_file, args.output_prefix, number_of_bins=args.number_of_bins,
                                            width_of_bins=args.width_of_bins, separator=args.separator,
                                            max_length=args.max_length, min_length=args.min_length,
                                            xlabel=args.xlabel, ylabel=args.ylabel, title=args.title,
                                            extensions=args.extensions, logbase=args.logbase)
"""
if (args.number_of_bins is not None) and (args.width_of_bins is not None):
    raise AttributeError("Options -w/--width_of_bins and -b/--number_of_bins mustn't be set simultaneously")

lengths = np.fromfile(args.input_file, sep=args.separator)

max_len = max(lengths)

if args.max_length is None:
    args.max_length = max_len

if (args.max_length != max_len) and (args.min_length != 1):
    filtered = []
    for entry in lengths:
        if args.min_length <= entry <= args.max_length:
            filtered.append(entry)
    else:
        filtered = lengths

figure = plt.figure(1, figsize=(6, 6))
subplot = plt.subplot(1, 1, 1)

if args.number_of_bins:
    bins = args.number_of_bins
elif args.width_of_bins:
    bins = np.arange(args.min_length, args.max_length, args.width_of_bins)
    #print bins
    #bins[0] += 1
    bins = np.append(bins, [args.max_length])
else:
    bins = 30

n, bins, patches = plt.hist(lengths, bins=bins)

bin_centers = (bins + ((bins[1] - bins[0])/2))[:-1]
#print bin_centers
#print len(n)
#print len(bin_centers)

plt.xlim(xmin=args.min_length, xmax=args.max_length)
if args.xlabel:
    plt.xlabel(args.xlabel)
if args.ylabel:
    plt.ylabel(args.ylabel)
if args.title:
    plt.title(args.title)

for ext in args.extensions:
    plt.savefig("%s.%s" % (args.output_prefix, ext))

subplot.set_yscale('log', basey=args.logbase)
#subplot.set_xscale('log', basex=args.logbase)
for ext in args.extensions:
    plt.savefig("%s.logscale.%s" % (args.output_prefix, ext))

# save histo values
np.savetxt("%s.histo" % args.output_prefix, zip(bin_centers, n), fmt="%i\t%i")
"""
