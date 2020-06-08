#!/usr/bin/env python3
__author__ = 'Sergei F. Kliver'

import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()

from RouToolPa.Routines import DrawingRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-f", "--filelist", action="store", dest="filelist", required=True, type=lambda s: s.split(","),
                    help="Comma-separated list of Id files(2-6)")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-l", "--set_labels", action="store", dest="set_labels", type=lambda s: s.split(","),
                    help="Comma-separated list of set labels")
parser.add_argument("-r", "--set_colors", action="store", dest="set_colors", type=lambda s: s.split(","),
                    help="Comma-separated list of set colors")

parser.add_argument("-e", "--extensions", action="store", dest="extensions", type=lambda x: x.split(","),
                    default=["png", "svg"],
                    help="Comma-separated list of extensions for histogram files")
parser.add_argument("-t", "--title", action="store", dest="title",
                    help="Title of histogram")

args = parser.parse_args()

DrawingRoutines.draw_venn_from_files(args.filelist, args.set_labels, args.output_prefix,
                                     extensions=args.extensions, title=args.title)


