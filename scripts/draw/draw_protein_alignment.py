#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

#from numpy import arange, int32, append, fromfile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
from RouToolPa.Collections.General import SynDict

from RouToolPa.Routines import DrawingRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_alignment", action="store", dest="alignment", required=True,
                    help="Input file with alignment")
parser.add_argument("-g", "--feature_gff", action="store", dest="features",
                    help="Gff file with features.Only features of first regions will be read")
parser.add_argument("-f", "--alignment_format", action="store", dest="alignment_format", default="fasta",
                    help="Format of alignment")
parser.add_argument("-e", "--extensions", action="store", dest="extensions", type=lambda x: x.split(","),
                    default=["png", "svg"],
                    help="Comma-separated list of extensions for picture")
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

"""

parser.add_argument("-s", "--separator", action="store", dest="separator", default="\n",
                    help="Separator between values in input file. Default - '\\n', i.e. one value per line")

parser.add_argument("-b", "--number_of_bins", action="store", dest="number_of_bins", type=int,
                    help="Number of bins in histogram. Incompatible with -w/--width_of_bins option. Default - 30")
parser.add_argument("-w", "--width_of_bins", action="store", dest="width_of_bins", type=float,
                    help="Width of bins in histogram. Incompatible with -b/--number_of_bins option. Not set by default")
parser.add_argument("-n", "--min_value", action="store", dest="min_length", type=float, default=0,
                    help="Minimum value to show. Default - 1")
parser.add_argument("-x", "--max_value", action="store", dest="max_length", type=float,
                    help="Maximum value to show. Default - length of longest sequence")


parser.add_argument("-l", "--xlabel", action="store", dest="xlabel",
                    help="X label")
parser.add_argument("-y", "--ylabel", action="store", dest="ylabel",
                    help="Y label")
parser.add_argument("-t", "--title", action="store", dest="title",
                    help="Title of histogram")
"""

args = parser.parse_args()

id_syn_dict = SynDict(filename=args.syn_file, key_index=args.key_column, value_index=args.value_column)

DrawingRoutines.draw_alignment_from_file(args.alignment, args.features, args.output_prefix,
                                         alignment_format=args.alignment_format, ext_list=args.extensions,
                                         label_fontsize=args.label_fontsize, left_offset=args.label_offset,
                                         figure_width=args.width, id_synonym_dict=id_syn_dict,
                                         id_replacement_mode=args.replacement_mode)
