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
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Directory to write output")
args = parser.parse_args()

NCBIRoutines.get_taxa_genomes_summary(args.taxa, args.email, args.output)
