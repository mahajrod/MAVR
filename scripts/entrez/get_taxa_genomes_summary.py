#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

from Bio import Entrez
from Routines.File import read_ids
from Routines import NCBIRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-t", "--taxa", action="store", dest="taxa", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of taxon names for query")
parser.add_argument("-e", "--email", action="store", dest="email", required=True,
                    help="Email used in Entrez queues")
parser.add_argument("-d", "--output_dir", action="store", dest="output_dir", default="./",
                    help="Directory to write output. Default: current directory")
parser.add_argument("-p", "--output_prefix", action="store", dest="output_prefix", default="",
                    help="Prefix of output files. Default: no prefix")
args = parser.parse_args()

NCBIRoutines.get_taxa_genomes_summary(args.taxa, args.email, args.output_dir, args.output_prefix)
