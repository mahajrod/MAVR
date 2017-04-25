#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

from Routines import DrawingRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-c", "--chromosome_gff", action="store", dest="chromosome_gff", required=True,
                    help="Input file with alignment")
parser.add_argument("-g", "--feature_gff", action="store", dest="feature_gff",
                    help="Gff file with features.")

parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    help="Prefix of output files")
parser.add_argument("-l", "--label_fontsize", action="store", dest="label_fontsize", default=13, type=int,
                    help="Label fontsize. Default: 13")
parser.add_argument("-a", "--label_offset", action="store", dest="label_offset", default=0.2, type=float,
                    help="Label offset. Increase if label is too long. Default: 0.2")
parser.add_argument("-w", "--width", action="store", dest="width", default=8, type=int,
                    help="Figure width(inches). Increase if label or alignment is too long. Default: 8")

parser.add_argument("-s", "--synonym_file", action="store", dest="syn_file",
                    help="File with id synonyms")
parser.add_argument("-k", "--key_column", action="store", dest="key_column", default=0, type=int,
                    help="Key column in synonym file(0-based). Default: 0")
parser.add_argument("-v", "--value_column", action="store", dest="value_column", default=1, type=int,
                    help="Value column in synonym file(1-based). Default: 1")
parser.add_argument("-r", "--replacement_mode", action="store", dest="replacement_mode", default="exact",
                    help="Replacement mode. Allowed: exact(default), partial")
parser.add_argument("-f", "--alignment_format", action="store", dest="alignment_format", default="fasta",
                    help="Format of alignment")
parser.add_argument("-e", "--extensions", action="store", dest="extensions", type=lambda x: x.split(","),
                    default=["png", "svg"],
                    help="Comma-separated list of extensions for picture")
parser.add_argument("-d", "--dpi", action="store", dest="dpi",
                    help="DPI of figure")

args = parser.parse_args()


DrawingRoutines.draw_chromosomes_with_features_simple(args.chromosome_gff, args.feature_gff, args.output_prefix,
                                                      figsize=(10, 10), dpi=args.dpi,
                                                      sense_feature_color="green", antisense_feature_color="red",
                                                      chromosome_color="black", label_fontsize=args.label_fontsize)
