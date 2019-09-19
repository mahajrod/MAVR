#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

from Pipelines.GenomeAnnotation import RepeatAnnotation

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input  fasta file with sequences")
parser.add_argument("-f", "--path_to_trf", action="store", dest="path_to_trf", default="trf",
                    help="Path to TRF binary")
parser.add_argument("-s", "--species", action="store", dest="species",
                    help="Species for RepeatMasker")
parser.add_argument("-l", "--library", action="store", dest="library",
                    help="Custom repeat library for RepeatMasker")
parser.add_argument("-d", "--output_dir", action="store", dest="output_dir", default="./",
                    help="Output directory. Default: ./")
parser.add_argument("-p", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Local(without /) prefix of output files")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int,
                    help="Number of threads")
parser.add_argument("-r", "--store_intermediate_file", action="store_true", dest="store",
                    help="Dont remove intermediate files")
args = parser.parse_args()

RepeatAnnotation.annotate_repeats(args.input, args.output_dir, args.output_prefix,
                                  repeatmasker=True, trf=True, windowmasker=True, threads=args.threads,
                                  trf_binary_path=args.path_to_trf,
                                  repeatmasker_soft_masking=True,
                                  repeatmasker_no_low_complexity=True,
                                  repeatmasker_custom_library=args.library,
                                  repeatmasker_species=args.species)
