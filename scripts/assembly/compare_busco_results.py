#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse
import matplotlib

matplotlib.use('Agg')
os.environ['MPLCONFIGDIR'] = '/tmp/'

import matplotlib.pyplot as plt
plt.ioff()
from RouToolPa.Tools.Annotation import BUSCO

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_busco_table_list", action="store", dest="input_busco_table_list", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of files with BUSCO results")
parser.add_argument("-l", "--label_list", action="store", dest="label_list",
                    type=lambda s: s.split(","),
                    help="Comma-separated list of assembly labels. Should have same length as list of "
                         "input files with assemblies. Default - not set, assemblies will be named like A1, A2, ../ ")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

args = parser.parse_args()

BUSCO.compare_multiple_genome_results(args.input_busco_table_list, args.output_prefix, label_list=args.label_list,
                                      black_scaffold_list=(), white_scaffold_list=())
