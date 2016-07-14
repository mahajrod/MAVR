#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

from Routines import MultipleAlignmentRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with alignment")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="File to write alignment of degenerate sites")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Alignment format. Default: fasta")
parser.add_argument("-g", "--genetic_code_table", action="store", dest="genetic_code_table", type=int,
                    default=1,
                    help="Genetic code table number")

args = parser.parse_args()

MultipleAlignmentRoutines.extract_degenerate_sites_from_codon_alignment_from_file(args.input, args.output,
                                                                                  genetic_code_table=args.genetic_code_table,
                                                                                  format=args.format)
