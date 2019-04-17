#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import MultipleAlignmentRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input alignment")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Output prefix")
parser.add_argument("-r", "--reference_sequence_id", action="store", dest="reference_sequence_id", required=True,
                    help="Reference sequence id")
parser.add_argument("-p", "--position list", action="store", dest="position_list", required=True,
                    type=lambda s: map(int, s.split(",")),
                    help="Comma separated list of positions in reference sequences")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of alignments")
parser.add_argument("-g", "--gap_symbol", action="store", dest="gap_symbol", default="-",
                    help="Gap symbol. Default - '-'")
parser.add_argument("-t", "--type", action="store", dest="type", default="nucleotide",
                    help="Alignment type. Allowed: nucleotide(default), codon, protein")
parser.add_argument("-l", "--flank_length", action="store", dest="flank_length", default=0, type=int,
                    help="Flank length. Default: 0, i.e no flanks will be included in the output file")

args = parser.parse_args()

MultipleAlignmentRoutines.get_specific_positions(args.input, args.reference_sequence_id, args.position_list,
                                                 args.output_prefix, format=args.format,
                                                 gap_symbol=args.gap_symbol, verbose=True,
                                                 alignment_type=args.type,
                                                 flank_length=args.flank_length)
