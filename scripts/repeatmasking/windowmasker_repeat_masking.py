#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Tools.BLAST import Windowmasker

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input_file", required=True,
                    help="Input file with sequences")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

args = parser.parse_args()


Windowmasker.masking(args.input_file, args.output_prefix, input_format="fasta", counts_format="obinary",
                     masking_format="interval", source="windowmasker", feature_type="repeat")


