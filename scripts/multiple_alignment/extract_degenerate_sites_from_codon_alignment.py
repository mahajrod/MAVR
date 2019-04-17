#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import MultipleAlignmentRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with alignment")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Alignment format. Default: fasta")
parser.add_argument("-g", "--genetic_code_table", action="store", dest="genetic_code_table", type=int,
                    default=1,
                    help="Genetic code table number")
parser.add_argument("-r", "--remove_Ns", action="store_true", dest="remove_Ns", default=False,
                    help="Remove codon columns with Ns. Default:False")

args = parser.parse_args()

MultipleAlignmentRoutines.extract_degenerate_sites_from_codon_alignment_from_file(args.input, args.output_prefix,
                                                                                  genetic_code_table=args.genetic_code_table,
                                                                                  format=args.format,
                                                                                  remove_codon_columns_with_Ns=args.remove_Ns)
