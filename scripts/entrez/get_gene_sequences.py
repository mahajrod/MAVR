#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

from Bio import Entrez
from Routines.File import read_ids
from Routines import NCBIRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-q", "--query", action="store", dest="query", required=True,
                    help="Query for gene database. it should contain all possible variants of gene names")
parser.add_argument("-e", "--email", action="store", dest="email", required=True,
                    help="Email used in Entrez queues")

args = parser.parse_args()

NCBIRoutines.get_gene_sequences(args.email, args.query)
