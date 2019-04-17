#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import EnsemblRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with protein sequences")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file")

args = parser.parse_args()

EnsemblRoutines.get_gene_transcript_protein_from_ensembl_pep_fasta(args.input, args.output)
