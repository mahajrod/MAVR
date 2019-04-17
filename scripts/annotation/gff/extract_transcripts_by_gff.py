#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Expression import Gffread

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_gff", action="store", dest="input_gff", required=True,
                    help="Input GFF file")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-g", "--genomic_fasta", action="store", dest="genomic_fasta", required=True,
                    help="Fasta with genome sequence")

args = parser.parse_args()

Gffread.extract_transcript_sequences(args.input_gff, args.genomic_fasta, args.output_prefix, coding_only=False)
