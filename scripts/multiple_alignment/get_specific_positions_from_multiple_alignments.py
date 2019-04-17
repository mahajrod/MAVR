#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import MultipleAlignmentRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_dir", action="store", dest="input_dir", required=True,
                    help="Input directory with alignments")
parser.add_argument("-o", "--output_dir", action="store", dest="output_dir", required=True,
                    help="Output directory")
parser.add_argument("-r", "--reference_sequence_id", action="store", dest="reference_sequence_id", required=True,
                    help="Reference sequence id")
parser.add_argument("-p", "--position_file", action="store", dest="position_file", required=True,
                    help="File with positions to extract")
parser.add_argument("-s", "--alignment_file_suffix", action="store", dest="alignment_file_suffix", default="",
                    help="Suffix of alignment files. Default: no suffix")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of alignments")
parser.add_argument("-g", "--gap_symbol", action="store", dest="gap_symbol", default="-",
                    help="Gap symbol. Default - '-'")
parser.add_argument("-t", "--type", action="store", dest="type", default="nucleotide",
                    help="Alignment type. Allowed: nucleotide(default), codon, protein")
parser.add_argument("-l", "--flank_length", action="store", dest="flank_length", default=0, type=int,
                    help="Flank length. Default: 0, i.e no flanks will be included in the output file")

args = parser.parse_args()

MultipleAlignmentRoutines.get_specific_positions_for_multiple_files(args.input_dir,
                                                                    args.position_file,
                                                                    args.reference_sequence_id,
                                                                    args.output_dir,
                                                                    alignment_file_suffix=args.alignment_file_suffix,
                                                                    format=args.format,
                                                                    gap_symbol=args.gap_symbol,
                                                                    verbose=True,
                                                                    alignment_type=args.type,
                                                                    flank_length=args.flank_length)
