#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse

from Bio import AlignIO
from Routines import MultipleAlignmentRoutines
from Routines.File import check_path, make_list_of_path_to_files, save_mkdir, split_filename


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", type=lambda x: x.split(","),
                    help="Comma-separated list of files or directory with files "
                         "containing alignments(one alignment per file)")
parser.add_argument("-n", "--max_gap_number", action="store", dest="max_gap_number", default=0, type=int,
                    help="Maximum number of gaps to retain column")
parser.add_argument("-o", "--output_directory", action="store", dest="output", type=check_path,
                    help="Output directory")
parser.add_argument("-g", "--gap_symbol", action="store", dest="gap_symbol", default="-",
                    help="Gap symbol")
parser.add_argument("-s", "--suffix", action="store", dest="suffix", default=".gaps_removed",
                    help="Suffix to use in output files. Default: '.gaps_removed'")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of alignment")


args = parser.parse_args()

save_mkdir(args.output)

for alignment_file in args.input:
    splited_filename = split_filename(alignment_file)
    output_filename = "%s%s%s%s" % (args.output, splited_filename[1], args.suffix, splited_filename[2])
    alignment = AlignIO.read(alignment_file, args.format)
    filtered_alignment = MultipleAlignmentRoutines.remove_columns_with_gaps(alignment, args.max_gap.number,
                                                                            gap_symbol=args.gap_symbol)
    AlignIO.write(filtered_alignment, output_filename, args.format)

