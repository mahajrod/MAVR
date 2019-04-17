#!/usr/bin/env python
__author__ = 'mahajrod'

import argparse
from RouToolPa.Tools.Bedtools import GenomeCov

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with per base coverage")
parser.add_argument("-s", "--sample_name", action="store", dest="sample_name", required=True,
                    help="Sample name for vcf file")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output vcf")
parser.add_argument("-r", "--reference", action="store", dest="reference", required=True,
                    help="Fasta file with reference")

args = parser.parse_args()

GenomeCov.convert_per_base_coverage_file_to_vcf(args.input, args.output, args.sample_name,
                                                args.reference,
                                                parsing_mode="parse", format="fasta")
