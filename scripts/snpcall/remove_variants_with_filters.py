#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse

from Tools.GATK import SelectVariants

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_vcf", action="store", dest="input_vcf", required=True,
                    help="Input vcf file")
parser.add_argument("-o", "--output_vcf", action="store", dest="output_vcf", required=True,
                    help="Output vcf file")
parser.add_argument("-r", "--reference", action="store", dest="reference", required=True,
                    help="Fasta with reference genome")

args = parser.parse_args()

SelectVariants.remove_entries_with_filters(args.reference, args.input_vcf, args.output_vcf)
