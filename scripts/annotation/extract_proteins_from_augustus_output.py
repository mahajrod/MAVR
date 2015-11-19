#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from Tools.Annotation import AUGUSTUS

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Gtf or gff file with augustus output")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="File to write proteins in fasta format")

args = parser.parse_args()

AUGUSTUS.extract_proteins_from_output(args.input, args.output)