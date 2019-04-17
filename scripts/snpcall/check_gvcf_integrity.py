#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import VCFRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_gvcf", action="store", dest="input_gvcf",
                    help="Input gvcf file",  required=True)
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-r", "--reference", action="store", dest="reference", required=True,
                    help="Fasta with reference genome")

args = parser.parse_args()

VCFRoutines.check_gvcf_integrity(args.input_gvcf,
                                 args.output_prefix,
                                 reference=args.reference,
                                 length_dict=None,
                                 parsing_mode="parse")
