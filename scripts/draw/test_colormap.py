#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
import matplotlib
from RouToolPa.Routines import DrawingRoutines
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()


parser = argparse.ArgumentParser()

parser.add_argument("-c", "--colormap", action="store", dest="colormap", required=True,
                    help="Name of matplotlib colormap to test")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    help="Prefix of output files")
parser.add_argument("-e", "--extensions", action="store", dest="extensions", type=lambda x: x.split(","),
                    default=["png"],
                    help="Comma-separated list of extensions for histogram files")
parser.add_argument("-n", "--color_number", action="store", dest="color_number", default=10,
                    type=int,
                    help="Number of colors")
parser.add_argument("-t", "--title", action="store", dest="title",
                    help="Title of histogram")

args = parser.parse_args()


DrawingRoutines.test_colormap(args.colormap, color_number=args.color_number,
                              output_prefix=args.output_prefix,
                              extension_list=args.extensions)
