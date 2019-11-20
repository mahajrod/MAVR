#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Tools.Population import PSMC

parser = argparse.ArgumentParser()

parser.add_argument("-s", "--sample_label_list", action="store", dest="sample_label_list", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of sample labels")
parser.add_argument("-p", "--psmc_list", action="store", dest="psmc_list", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of PSMC files corresponding to sample labels")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-t", "--generation_time", action="store", dest="generation_time",
                    required=True, type=int,
                    help="Generation time in years(int)")
parser.add_argument("-a", "--mutation_rate", action="store", dest="mutation_rate", required=True, type=float,
                    help="Absolute mutation rate per year per nucleotide")
parser.add_argument("-m", "--min_threshold", action="store", dest="min_threshold", type=float,
                    help="Minimum threshold for generations to plot. Use 0 for auto. Default: 10000")
parser.add_argument("-x", "--max_threshold", action="store", dest="max_threshold", type=float,
                    help="Maximum threshold for generations to plot. Default: auto")

parser.add_argument("-g", "--plot_grid", action="store_true", dest="plot_grid",
                    help="Plot grid. Default: False")

parser.add_argument("-l", "--psmc_dir", action="store", dest="psmc_dir", default="",
                    help="Path to directory with psmc_plot.pl script")

args = parser.parse_args()

PSMC.path = args.psmc_dir
PSMC.psmc_plot(args.sample_label_list, args.psmc_list, args.generation_time,
               args.mutation_rate, args.output_prefix, plot_grid=args.plot_grid,
               min_generations=args.min_threshold, max_generations=args.max_threshold)
