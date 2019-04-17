#!/usr/bin/env python
__author__ = 'mahajrod'

import argparse
from RouToolPa.Tools.GATK import FastaAlternateReferenceMaker

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with database")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file")

args = parser.parse_args()

FastaAlternateReferenceMaker.convert_picard_dict_to_syn_file(args.input, args.output)
