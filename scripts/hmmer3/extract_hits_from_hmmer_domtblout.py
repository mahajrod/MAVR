#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.HMMER import HMMER3

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_hmm", action="store", dest="input", required=True,
                    help="Input hmmer3 domtblout file")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Output prefix")

args = parser.parse_args()

output_pfam_annotated_dom_ids = "%s.dom_ids" % args.output_prefix
output_pfam_annotated_dom_names = "%s.dom_names" % args.output_prefix

HMMER3.extract_dom_ids_hits_from_domtblout(args.input, output_pfam_annotated_dom_ids)
HMMER3.extract_dom_names_hits_from_domtblout(args.input, output_pfam_annotated_dom_names)
