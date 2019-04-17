#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import SequenceRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-f", "--forward", action="store", dest="forward", required=True,
                    help="Comma separated list of files/directories with forward sequences")
parser.add_argument("-r", "--reverse", action="store", dest="reverse", required=True,
                    help="Comma separated list of files/directories with reverse sequences")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

parser.add_argument("-a", "--forward_suffix", action="store", dest="forward_suffix", required=True,
                    help="Suffix of sequence ids from forward files")
parser.add_argument("-b", "--reverse_suffix", action="store", dest="reverse_suffix", required=True,
                    help="Suffix of sequence ids from reverse files")

parser.add_argument("-m", "--format", action="store", dest="format", default="fasta",
                    help="Format of input and output file. Allowed formats genbank, fasta(default)")
parser.add_argument("-p", "--parsing_mode", action="store", dest="parsing_mode", default="index_db",
                    help="Parsing mode for input sequence file. "
                         "Possible variants: 'index_db'(default), 'index', 'parse'")

args = parser.parse_args()


"""
example of usage

~/Soft/MAVR/scripts/sequence/check_pairing.py -p parse \
                                              -a ".F" \
                                              -b ".R" \
                                              -o GSS_BOH_BAC_end \
                                              -f GSS_BOH_BAC_end.forward.fa \
                                              -r GSS_BOH_BAC_end.reverse.fa

"""
SequenceRoutines.check_pairing_from_file(args.forward, args.reverse, args.output_prefix, args.forward_suffix,
                                         args.reverse_suffix, parsing_mode=args.parsing_mode, format=args.format,
                                         forward_index_file="forward_tmp.idx", reverse_index_file="reverse_tmp.idx",
                                         retain_index=False, output_file_extension=args.format)


