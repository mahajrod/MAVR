#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import DrawingRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-c", "--chromosome_gff", action="store", dest="chromosome_gff", required=True,
                    help="Input file with alignment")
parser.add_argument("-g", "--feature_gff", action="store", dest="feature_gff",
                    help="Gff file with features.")

parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    help="Prefix of output files")
parser.add_argument("-n", "--id_field", action="store", dest="id_field", default="ID",
                    help="Id field in feature gff file. Default: ID")
parser.add_argument("-r", "--chromosome_name_field", action="store", dest="chr_field", default=None,
                    help="Chromosome name field in chromosome gff file. Default: use chromosome id")
parser.add_argument("-l", "--label_fontsize", action="store", dest="label_fontsize", default=13, type=int,
                    help="Label fontsize. Default: 13")
parser.add_argument("-a", "--label_offset", action="store", dest="label_offset", default=0.2, type=float,
                    help="Label offset. Increase if label is too long. Default: 0.2")
parser.add_argument("-w", "--width", action="store", dest="width", default=8, type=int,
                    help="Figure width(inches). Increase if label or alignment is too long. Default: 8")

parser.add_argument("-p", "--deseq2_pwc_file", action="store", dest="deseq2_pwc_file",
                    help="File results of DESeq2 pairwise comparison")

parser.add_argument("-e", "--extensions", action="store", dest="extensions", type=lambda x: x.split(","),
                    default=["png", "svg"],
                    help="Comma-separated list of extensions for picture")
parser.add_argument("-d", "--dpi", action="store", dest="dpi", type=int,
                    help="DPI of figure")

args = parser.parse_args()

DrawingRoutines.draw_chromosomes_with_features_simple(args.chromosome_gff, args.feature_gff, args.output_prefix,
                                                      figsize=(10, 10), dpi=args.dpi, ext_list=args.extensions,
                                                      id_field_in_gff=args.id_field,
                                                      chromosome_name_field=args.chr_field,
                                                      sense_feature_color="green", antisense_feature_color="red",
                                                      chromosome_color="black", label_fontsize=args.label_fontsize,
                                                      deseq2_pwc_file=args.deseq2_pwc_file,
                                                      upregulated_color="grey", downregulated_color="lightgrey",
                                                      #upregulated_color="green", downregulated_color="red",
                                                      absent_expression_data_color="black",
                                                      coloring_mode="strand" if args.deseq2_pwc_file is None else "expression"
                                                      )
