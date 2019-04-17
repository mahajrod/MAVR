#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import MatplotlibRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    type=lambda s : s.split(","),
                    help="Comma separated list of two input files with data")
parser.add_argument("-d", "--index", action="store", dest="index",
                    type=lambda s: map(int, s.split(",")),
                    help="Zero based indexes of file columns to use. Default: all ")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    help="Prefix of output files")
parser.add_argument("-s", "--separator", action="store", dest="separator", default=None,
                    help="Separator between values in input file. Default - any space symbol")
parser.add_argument("-l", "--xlabel", action="store", dest="xlabel", type=lambda s: s.split(","),
                    help="Comma-separated list of X labels")
parser.add_argument("-y", "--ylabel", action="store", dest="ylabel", type=lambda s: s.split(","),
                    help="Comma-separated list of Y labels")
parser.add_argument("-b", "--number_of_bins", action="store", dest="number_of_bins",
                    type=lambda s: map(int, s.split(",")),
                    help="Comma-separated list of bin numbers in histograms. "
                         "Incompatible with -w/--width_of_bins option. Default - 30")

parser.add_argument("-w", "--width_of_bins", action="store", dest="width_of_bins",
                    type=lambda s: map(float, s.split(",")),
                    help="Comma-separated list of bin widths in histograms. "
                         "Incompatible with -b/--number_of_bins option. Not set by default")
parser.add_argument("-n", "--min_value", action="store", dest="min_value", type=lambda s: map(float, s.split(",")),
                    default=0,
                    help="Comma-separated list of minimum value to show. Default - 1")
parser.add_argument("-x", "--max_value", action="store", dest="max_value", type=lambda s: map(float, s.split(",")),
                    help="Comma-separated list of maximum value to show. Default - length of longest sequence")
parser.add_argument("-g", "--logbase", action="store", dest="logbase", type=int, default=10,
                    help="Logbase to use for log-scaled histograms")
parser.add_argument("-e", "--extensions", action="store", dest="extensions", type=lambda x: x.split(","),
                    default=["png"],
                    help="Comma-separated list of extensions for histogram files. Default: png only")

parser.add_argument("-t", "--title_list", action="store", dest="title_list", type=lambda s: s.split(","),
                    help="Comma-separated ist of two title for histograms")
parser.add_argument("-v", "--share_x_axis", action="store_true", dest="share_x_axis", default=False,
                    help="Share X axis in all histograms. Default: False")
parser.add_argument("-z", "--share_y_axis", action="store_true", dest="share_y_axis", default=False,
                    help="Share Y axis in all histograms. Default: False")

args = parser.parse_args()

if args.index is None:
    args.index = [None for i in range(0, len(args.input))]
if args.max_value is None:
    args.max_value = [None for i in range(0, len(args.input))]

MatplotlibRoutines.draw_tetra_histogram_with_two_logscaled_from_file(args.input, args.index, args.output_prefix,
                                                                     figsize=(10, 10),
                                                                     number_of_bins_list=args.number_of_bins,
                                                                     width_of_bins_list=args.width_of_bins,
                                                                     max_threshold_list=args.max_value,
                                                                     min_threshold_list=args.min_value,
                                                                     xlabel=args.xlabel, ylabel=args.ylabel,
                                                                     title_list=args.title_list,
                                                                     logbase=args.logbase, label_list=None,
                                                                     extensions=args.extensions, suptitle=None,
                                                                     separator=args.separator,
                                                                     share_y_axis=args.share_y_axis,
                                                                     share_x_axis=args.share_x_axis)
"""
Example:
~/Dropbox/MAVR/scripts/draw/draw_tetra_histogram_with_two_logscaled.py -i kirill.dn.ds.w.tab,solenodon.raw_alns.all.tab -d 3,3 -o dnds.ratio.log  -l 'dN/dS' -y "Number of genes" -w 20 -n 0 -x 999 -t "11 species,4 species"
"""
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
