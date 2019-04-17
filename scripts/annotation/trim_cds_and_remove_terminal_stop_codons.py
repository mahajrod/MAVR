#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import SequenceRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input fasta with raw cDS")
parser.add_argument("-s", "--stop_codons", action="store", dest="stop_codons_list", default=["TGA", "TAA", "TAG"],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of stop codons. Can be set using any case and both RNA and DNA alphabet."
                         "Default: TGA, TAA, TAG")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file to write trimmed CDS")

args = parser.parse_args()
print("Using %s as stop codons" % ",".join(args.stop_codons_list))
SequenceRoutines.trim_cds_and_remove_terminal_stop_codons(args.input, args.output, args.stop_codons_list)
