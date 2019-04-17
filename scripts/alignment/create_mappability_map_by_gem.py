#!/usr/bin/env python
__author__ = 'mahajrod'

import argparse
from RouToolPa.Tools.Alignment import GEM

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--gem_index", action="store", dest="gem_index",
                    help="Gem index")
parser.add_argument("-r", "--reference", action="store", dest="reference",
                    help="Fasta file with reference")
parser.add_argument("-k", "--kmer_list", action="store", dest="kmer_list", default=[50, 75, 100, 125, 150],
                    type=lambda s: map(int, GEM.split_string_by_comma(s)),
                    help="Comma-separated list of kmers. Default: 50, 75, 100, 125, 150")

parser.add_argument("-g", "--gem_dir", action="store", dest="gem_dir",
                    default="", help="Path to directory with GEM binaries")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=4,
                    help="Number of threads")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

args = parser.parse_args()

GEM.threads = args.threads
GEM.path = args.gem_dir
GEM.create_map_of_mappability(args.kmer_list, args.output_prefix, gem_index=args.gem_index,
                              reference=args.reference, convert_to_wig=True)
