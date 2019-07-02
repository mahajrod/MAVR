#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os
import argparse
from Bio import SeqIO

os.environ['MPLCONFIGDIR'] = '/tmp/'
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()

from RouToolPa.Routines.File import make_list_of_path_to_files
from RouToolPa.Tools.Kmers import Jellyfish
from RouToolPa.Routines.Sequence import rev_com_generator

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", type=lambda s: s.split(","),
                    help="Comma-separated list of fasta or fastq files or directories containing them.")
parser.add_argument("-m", "--kmer_length", action="store", dest="kmer_length", type=int, default=23,
                    help="Length of kmers")
parser.add_argument("-s", "--hash_size", action="store", dest="hash_size", type=int, default=1000000,
                    help="Size of hash. Estimation of hash size: for short reads S=(G + k*n)/0.8, "
                    "G - genome size, k - kmer length, n - number of reads, for assembled sequences "
                    "S=Sum(L)")
parser.add_argument("-a", "--base_prefix", action="store", dest="base_prefix", default="jellyfish_db",
                    help="Name of kmer database. Default: jellyfish_db")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads")
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
parser.add_argument("-d", "--draw_distribution", action="store_true", dest="draw_distribution",
                    help="Draw distribution of kmers")
parser.add_argument("-j", "--jellyfish_path", action="store", dest="jellyfish_path",
                    help="Path to jellyfish")
parser.add_argument("-n", "--dont_extract_kmer_list", action="store_true", dest="dont_extract_kmer_list",
                    help="Don't extract kmer list")

args = parser.parse_args()

args.input = make_list_of_path_to_files(args.input)
if args.count_both_strands and args.add_rev_com:
    raise ValueError("Options -b/--count_both_strands and -r/--add_reverse_complement are not compatible")

if args.add_rev_com:
    file_with_rev_com = args.base_prefix + "_with_rev_com.fasta"
    record_dict = SeqIO.index_db("temp_index.idx", args.input, format="fasta")
    SeqIO.write(rev_com_generator(record_dict, yield_original_record=True), file_with_rev_com, "fasta")
    args.base_prefix += "_with_rev_com"

base_file = "%s_%i_mer.jf" % (args.base_prefix, args.kmer_length)
kmer_table_file = "%s_%i_mer.counts" % (args.base_prefix, args.kmer_length)
kmer_file = "%s_%i_mer.kmer" % (args.base_prefix, args.kmer_length)

Jellyfish.threads = args.threads
Jellyfish.path = args.jellyfish_path if args.jellyfish_path else ""
Jellyfish.count(args.input if not args.add_rev_com else file_with_rev_com, base_file,
                kmer_length=args.kmer_length, hash_size=args.hash_size,
                count_both_strands=args.count_both_strands)
if not args.dont_extract_kmer_list:
    Jellyfish.dump(base_file, kmer_table_file)
    sed_string = 'sed -e "s/\t.*//" %s > %s' % (kmer_table_file, kmer_file)
    os.system(sed_string)

if args.draw_distribution:
    histo_file = "%s_%i_mer.histo" % (args.base_prefix, args.kmer_length)
    picture_prefix = "%s_%i_mer_histogram" % (args.base_prefix, args.kmer_length)
    Jellyfish.histo(base_file, histo_file, upper_count=10000000)

    counts = []
    bins = []

    with open(histo_file, "r") as histo_fd:
        for line in histo_fd:
            entry = line.strip().split()
            counts.append(entry[1])
            bins.append(entry[0])

    figure = plt.figure(1, figsize=(8, 8), dpi=300)
    subplot = plt.subplot(1, 1, 1)
    plt.suptitle("Distribution of %i-mers" % args.kmer_length,
                 fontweight='bold')
    plt.plot(bins, counts)

    plt.xlabel("Multiplicity")
    plt.ylabel("Number of distinct kmers")
    subplot.set_yscale('log', basey=10)
    subplot.set_xscale('log', basex=10)
    for extension in ["png", "svg"]:
        plt.savefig("%s.%s" % (picture_prefix, extension))

if args.add_rev_com:
    os.remove("temp_index.idx")
