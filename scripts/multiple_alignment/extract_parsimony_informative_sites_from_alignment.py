#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import MultipleAlignmentRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with alignment")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Alignment format. Default: fasta")
parser.add_argument("-r", "--remove_Ns", action="store_true", dest="remove_Ns", default=False,
                    help="Remove columns with Ns. Default:False")
parser.add_argument("-a", "--remove_columns_with_ambigous_nucleotides", action="store_true",
                    dest="remove_columns_with_ambigous_nucleotides", default=False,
                    help="Remove columns with a. Default:False")

args = parser.parse_args()

MultipleAlignmentRoutines.extract_parsimony_informative_sites_from_alignment_from_file(args.input,
                                                                                       args.output,
                                                                                       format=args.format,
                                                                                       remove_columns_with_Ns=args.remove_Ns,
                                                                                       remove_columns_with_ambigous_nucleotides=args.remove_columns_with_ambigous_nucleotides)
