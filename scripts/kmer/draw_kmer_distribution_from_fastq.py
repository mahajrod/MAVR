#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os
import argparse

#import numpy as np
#import matplotlib
#matplotlib.use('Agg')
#os.environ['MPLCONFIGDIR'] = '/tmp/'
#import matplotlib.pyplot as plt

from Bio import SeqIO

#from RouToolPa.Routines import MatplotlibRoutines
from RouToolPa.Routines.Sequence import rev_com_generator
from RouToolPa.Routines.File import make_list_of_path_to_files

from RouToolPa.Tools.Kmers import Jellyfish


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", type=lambda s: s.split(","), required=True,
                    help="Comma-separated list of fasta or fastq files or directories containing them.")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Output prefix")
parser.add_argument("-e", "--output_formats", action="store", dest="output_formats", type=lambda s: s.split(","),
                    default=["svg", "png"],
                    help="Comma-separated list of formats (supported by matlotlib) "
                         "of output figure.Default: svg,eps,pdf,png,jpg")
parser.add_argument("-l", "--logbase", action="store", dest="logbase", type=int, default=10,
                    help="Base of logarithm. Default - 10")
parser.add_argument("-m", "--kmer_length", action="store", dest="kmer_length", type=int, default=23,
                    help="Length of kmers. Default - 23")
parser.add_argument("-s", "--hash_size", action="store", dest="hash_size", type=str, default="1000000",
                    help="Size of hash. Estimation of hash size: for short reads S=(G + k*n)/0.8, "
                    "G - genome size, k - kmer length, n - number of reads, for assembled sequences "
                    "S=Sum(L). Default - 1000000")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads. Default - 1")
parser.add_argument("-b", "--count_both_strands", action="store_true", dest="count_both_strands",
                    help="Count kmers in both strands. NOTICE: only mer or its reverse-complement, whichever "
                         "comes first lexicographically, is stored and the count value is the number of "
                         "occurrences of both. So this option is not suitable for generating sets of forward "
                         "and reverse-complement kmers. For this case use -r/--add_reverse_complement option. "
                         "Not compatible with -r/--add_reverse_complement option.")
parser.add_argument("-r", "--add_reverse_complement", action="store_true", dest="add_rev_com",
                    help="Add reverse-complement sequences before counting kmers. "
                         "Works only for fasta sequences. "
                         "Not compatible with -b/--count_both_strands option")
parser.add_argument("-j", "--jellyfish_path", action="store", dest="jellyfish_path",
                    help="Path to jellyfish")
parser.add_argument("-w", "--low_limit", action="store", dest="low_limit", type=int, default=5,
                    help="Low limit of histogram without logscale. Default - 5")
parser.add_argument("-g", "-high_limit", action="store", dest="high_limit", type=int, default=100,
                    help="High limit of histogram without logscale. Default - 100")
#parser.add_argument("-d", "--draw_peaks_and_gaps", action="store_true", dest="draw_peaks_and_gaps",
#                    help="Draw peaks and gaps")

args = parser.parse_args()

args.input = make_list_of_path_to_files(args.input)
if args.count_both_strands and args.add_rev_com:
    raise ValueError("Options -b/--count_both_strands and -r/--add_reverse_complement are not compatible")

if args.add_rev_com:
    file_with_rev_com = args.output_prefix + "_with_rev_com.fasta"
    record_dict = SeqIO.index_db("temp_index.idx", args.input, format="fasta")
    SeqIO.write(rev_com_generator(record_dict, yield_original_record=True), file_with_rev_com, "fasta")
    args.output_prefix += "_with_rev_com"
    os.remove("temp_index.idx")

base_file = "%s_%i_mer.jf" % (args.output_prefix, args.kmer_length)
kmer_table_file = "%s_%i_mer.counts" % (args.output_prefix, args.kmer_length)
kmer_file = "%s_%i_mer.kmer" % (args.output_prefix, args.kmer_length)

histo_file = "%s_%i_mer.histo" % (args.output_prefix, args.kmer_length)
picture_prefix = "%s_%i_mer_histogram" % (args.output_prefix, args.kmer_length)

Jellyfish.threads = args.threads
Jellyfish.timelog = "%s_%i_mer.jellyfish.time.log" % (args.output_prefix, args.kmer_length)
Jellyfish.path = args.jellyfish_path if args.jellyfish_path else ""
Jellyfish.count(args.input if not args.add_rev_com else file_with_rev_com, base_file,
                kmer_length=args.kmer_length, hash_size=args.hash_size,
                count_both_strands=args.count_both_strands)
Jellyfish.histo(base_file, histo_file, upper_count=10000000)
Jellyfish.draw_kmer_distribution(histo_file, args.kmer_length, picture_prefix, output_formats=args.output_formats,
                                 logbase=args.logbase, non_log_low_limit=args.low_limit,
                                 non_log_high_limit=args.high_limit) #, draw_peaks_and_gaps=args.draw_peaks_and_gaps)
"""
bins, counts = np.loadtxt(args.input, unpack=True)

#counts = []
#bins = []
#
#with open(histo_file, "r") as histo_fd:
#    for line in histo_fd:
#        entry = line.strip().split()
#        counts.append(entry[1])
#        bins.append(entry[0])

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
    plt.savefig("%s.logscale.%s" % (picture_prefix, extension))

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
    plt.savefig("%s.no_logscale.%s" % (picture_prefix, extension))

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
    plt.savefig("%s.%s" % (picture_prefix, extension))
"""
