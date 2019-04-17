#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Alignment import LongRanger



parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="File with input reference sequences")
parser.add_argument("-o", "--output", action="store", dest="output",
                    help="File to write prepared reference")
parser.add_argument("-c", "--coord_file", action="store", dest="coord_file",
                    help="File to write coordinates of sequences in merged record")
parser.add_argument("-l", "--longranger_dir", action="store", dest="longranger_dir", default="",
                    help="Directory with longranger binary")
parser.add_argument("-e", "--renamed_scaffolds_syn_file", action="store", dest="renamed_scaffolds_syn_file",
                    required=True,
                    help="File with synonyms to renamed scaffolds")
parser.add_argument("-r", "--symbols_to_replace_list", action="store", dest="symbols_to_replace_list", default=[":"],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of symbols to replace in scaffold ids. Default: :")
args = parser.parse_args()


LongRanger.path = args.longranger_dir
LongRanger.prepare_reference(args.input, args.output, args.renamed_scaffolds_syn_file,
                             max_scaffold_length=527000000, max_scaffold_number=500,
                             polyN_len=500, coord_file=args.coord_file, symbols_to_replace_list=args.symbols_to_replace_list)
