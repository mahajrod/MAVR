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

parser.add_argument("-i", "--input_fasta_list", action="store", dest="input_fasta_list", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of fasta files with different assemblies")
parser.add_argument("-l", "--label_list", action="store", dest="label_list",
                    type=lambda s: s.split(","),
                    help="Comma-separated list of assembly labels. Should have same length as list of "
                         "input files with assemblies. Default - not set, assemblies will be named like A1, A2, ../ ")

parser.add_argument("-d", "--busco_db", action="store", dest="busco_db", required=True,
                    help="BUSCO database")
parser.add_argument("-s", "--species", action="store", dest="species", default="human",
                    help="Species for AUGUSTUS. Default: human")
parser.add_argument("-o", "--output_dir", action="store", dest="output_dir", default="./",
                    help="Output directory. Default: ./")

parser.add_argument("-m", "--mode", action="store", dest="mode", default="genome",
                    help="BUSCO mode. Allowed: genome(default), transcriptome, proteins")

parser.add_argument("-u", "--busco_dir", action="store", dest="busco_dir",
                    help="Path to directory with BUSCO script")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=4,
                    help="Number of threads to use")
args = parser.parse_args()

BUSCO.path = args.busco_dir
BUSCO.threads = args.threads
BUSCO.assess_multiple_genomes(args.input_fasta_list, args.busco_db, args.species,
                              label_list=args.label_list,
                              output_dir=args.output_dir, mode=args.mode)
