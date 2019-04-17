#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Tools.Phylogenetics import Gblocks


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input_dir", required=True,
                    help="Input directory with alignments")
parser.add_argument("-o", "--output_dir", action="store", dest="output_dir", required=True,
                    help="Output directory")
parser.add_argument("-p", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Output prefix")
parser.add_argument("-y", "--input_type", action="store", dest="input_type", default="codon",
                    help="Type of input sequences. Allowed: dna, codon(default), protein")
parser.add_argument("-m", "--min_seq_number_for_conserved_position", action="store",
                    dest="min_seq_number_for_conserved_position", type=int,
                    help="Minimum number of sequences to treat position conservative")
parser.add_argument("-f", "--min_seq_number_for_flank_position", action="store",
                    dest="min_seq_number_for_flank_position", type=int,
                    help="Minimum number of sequences for flank positions")
parser.add_argument("-x", "--max_pos_number_for_noncons_contig_pos", action="store",
                    dest="max_pos_number_for_noncons_contig_pos", type=int,
                    help="Maximum number of contiguous nonconserved positions")
parser.add_argument("-n", "--min_block_len", action="store",dest="min_block_len", type=int,
                    help="Minimum block length")
parser.add_argument("-t", "--threads", action="store",dest="threads", type=int, default=1,
                    help="Thread number")

args = parser.parse_args()

Gblocks.threads = args.threads

Gblocks.parallel_run(args.input_dir, args.output_dir, args.output_prefix,
                     input_type=args.input_type,
                     min_seq_number_for_conserved_position=args.min_seq_number_for_conserved_position,
                     min_seq_number_for_flank_position=args.min_seq_number_for_flank_position,
                     max_pos_number_for_noncons_contig_pos=args.max_pos_number_for_noncons_contig_pos,
                     min_block_len=args.min_block_len,
                     allow_gaps="half",
                     save_postscript=True,
                     output_type="htm",
                     )
