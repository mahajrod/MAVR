#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import MultipleAlignmentRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="File with alignment")
parser.add_argument("-c", "--coordinates", action="store", dest="coordinates", required=True,
                    help="File with coordinates of gene alignments")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="File to write alignment prepared for Codeml")

args = parser.parse_args()

MultipleAlignmentRoutines.prepare_multigene_alignment_for_codeml(args.input, args.coordinates, args.output,
                                                                 format="fasta")
