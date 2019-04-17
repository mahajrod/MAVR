#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.GATK import CombineVariants


parser = argparse.ArgumentParser()

parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output vcf file")
parser.add_argument("-r", "--reference", action="store", dest="reference", required=True,
                    help="Fasta with reference genome")
parser.add_argument("-g", "--gatk_directory", action="store", dest="gatk_dir", default="",
                    help="Directory with GATK jar")
parser.add_argument("-i", "--input_vcf_list", action="store", dest="input_vcf_list", type=lambda s: s.split(","),
                    help="Comma-separated list of vcf files to combine",  required=True,)

args = parser.parse_args()


CombineVariants.combine_from_same_source(args.gatk_dir, args.reference, args.input_vcf_list, args.output)
