#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

from Tools.Population import PSMC

parser = argparse.ArgumentParser()

parser.add_argument("-s", "--sample_label_list", action="store", dest="sample_label_list", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of sample labels")
parser.add_argument("-p", "--psmc_list", action="store", dest="psmc_list", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of PSMC files corresponding to sample labels")
parser.add_argument("-t", "--generation_time", action="store", dest="generation_time",
                    required=True, type=int,
                    help="Generation time in years(int)")
parser.add_argument("-a", "--mutation_rate", action="store", dest="mutation_rate", required=True, type=float,
                    help="Absolute mutation rate per year per nucleotide")
parser.add_argument("-g", "--plot_grid", action="store_true", dest="plot_grid",
                    help="Plot grid. Default: False")
args = parser.parse_args()

PSMC.psmc_plot(args.sample_label_list, args.psmc_list, args.generation_time,
               args.mutation_rate, plot_grid=args.plot_grid)
