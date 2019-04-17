#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.MultipleAlignment import PRANK
from RouToolPa.Routines import FileRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    type=lambda x: FileRoutines.make_list_of_path_to_files(x.split(",")),
                    help="Comma-separated list of files or directory with files containing sequences to be aligned")
parser.add_argument("-p", "--processes", action="store", dest="processes", type=int, default=1,
                    help="Number of simultaneously running alignments")
parser.add_argument("-o", "--output_directory", action="store", dest="output", type=FileRoutines.check_path,
                    required=True,
                    help="Output directory")
parser.add_argument("-t", "--tree_file", action="store", dest="tree_file",
                    help="File with tree")
parser.add_argument("-s", "--suffix", action="store", dest="suffix", default="",
                    help="Suffix of basename of output_files")
parser.add_argument("-e", "--skip_insertions", action="store_true", dest="skip_insertions",
                    help="Skip insertions")

args = parser.parse_args()

PRANK.threads = args.processes
PRANK.parallel_align(args.input, args.output, output_suffix=args.suffix, tree_file=args.tree_file, output_format=None,
                     show_xml=None,
                     show_tree=None, show_ancestral_sequences=None, show_evolutionary_events=None,
                     showall=None, compute_posterior_support=None, njtree=None, skip_insertions=args.skip_insertions,
                     codon_alignment=None, translated_alignment=None)

