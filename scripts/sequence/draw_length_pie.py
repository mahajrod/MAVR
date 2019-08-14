#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse


from RouToolPa.Routines import SequenceRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input",
                    help="Input file with scaffold ids and its lengths")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    help="Prefix of output files")
parser.add_argument("-e", "--extensions", action="store", dest="extensions", type=lambda x: x.split(","),
                    default=["png"],
                    help="Comma-separated list of extensions for histogram files")


args = parser.parse_args()

SequenceRoutines.draw_length_pie(args.input, args. output_prefix, thresholds=(1000, 10000, 100000, 1000000, 10000000),
                                 labels=("1-10 kbp", "10-100 kbp", "0.1-1 Mbp", "1-10 Mbp", "10+ Mbp"),
                                 colors=("red", "orange", "yellow", "blue", "green"),
                                 extension_list=args.extensions)

