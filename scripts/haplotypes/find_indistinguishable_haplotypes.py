#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import HaplotypeRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--seq_file", action="store", dest="seq_file", required=True,
                    help="Input file with sequences")
parser.add_argument("-s", "--haplotype_file", action="store", dest="haplotype_file",
                    help="File with haplotype synonyms to sequence ids. If not set sequence ids will be used.")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-t", "--threads", action="store", dest="threads", default=1, type=int,
                    help="Thread number to use for clustering. Default: 1")

parser.add_argument("-c", "--cdhit_dir", action="store", dest="cdhit_dir", default=None,
                    help="Path to CDHIT directory")

args = parser.parse_args()

HaplotypeRoutines.find_indistinguishable_haplotypes(args.seq_file, args.haplotype_file, args.output_prefix,
                                                    threads=args.threads, cdhit_dir=args.cdhit_dir)
