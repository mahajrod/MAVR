#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import HaplotypeRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-a", "--alignment", action="store", dest="aln_file", required=True,
                    help="Input file with alignment in fasta format")
parser.add_argument("-y", "--hap_file", action="store", dest="hap_file", required=True,
                    help="Fam file with haplotypes")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file")


args = parser.parse_args()

HaplotypeRoutines.prepare_template_for_popart(args.aln_file, args.hap_file, args.output)
