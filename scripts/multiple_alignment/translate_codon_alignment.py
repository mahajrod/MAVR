#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import MultipleAlignmentRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-c", "--codon_alignment", action="store", dest="codon_alignment", required=True,
                    help="Input file with codon alignment")
parser.add_argument("-p", "--protein_alignment", action="store", dest="protein_alignment", required=True,
                    help="Output file with protein alignment")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of alignments. Default: fasta")
parser.add_argument("-g", "--gap_symbol", action="store", dest="gap_symbol", default="-",
                    help="Gap symbol in alignment. Default: '-'")
parser.add_argument("-t", "--genetic_code", action="store", dest="genetic_code", default=1, type=int,
                    help="Genetic code to use(NCBI tables) . Default: 1(standart)")

args = parser.parse_args()

MultipleAlignmentRoutines.translate_codon_alignment(args.codon_alignment, args.protein_alignment, format=args.format,
                                                    gap_symbol=args.gap_symbol, table=args.genetic_code)
