#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Tools.Population import PLINK

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--plink_binary_data_prefix", action="store", dest="plink_binary_data_prefix", required=True,
                    help="Prefix of PLINK binary data. Files with.fam, .bed and .bim are expected")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file with sequences")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of output file. Default: fasta")
parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                    help="Verbose")

args = parser.parse_args()

PLINK.extract_sequences_from_plink_binary_snp_data(args.plink_binary_data_prefix, args.output, verbose=args.verbose,
                                                   output_format=args.format)
