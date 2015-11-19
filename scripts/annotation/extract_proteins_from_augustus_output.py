#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from Tools.Annotation import AUGUSTUS

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Gtf or gff file with augustus output")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="File to write proteins in fasta format")
parser.add_argument("-d", "--id_prefix", action="store", dest="id_prefix", default="p.",
                    help="Prefix to use for protein ids")
args = parser.parse_args()

AUGUSTUS.extract_proteins_from_output(args.input, args.output, id_prefix=args.id_prefix)