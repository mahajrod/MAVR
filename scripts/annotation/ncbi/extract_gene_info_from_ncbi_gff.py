#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse

from RouToolPa.Routines import NCBIRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-g", "--input_gff", action="store", dest="input_gff", required=True,
                    help="Input gff with NCBI annotations")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file")
parser.add_argument("-f", "--feature_type", action="store", dest="feature_type",
                    type=lambda s: s.split(","), default=["gene"],
                    help="Feature type to use")


args = parser.parse_args()

NCBIRoutines.extract_gene_info_from_ncbi_gff(args.input_gff, args.output, feature_type_list=args.feature_type)
