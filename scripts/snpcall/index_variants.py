#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.GATK4 import ValidateVariants4


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_vcf", action="store", dest="input_vcf", required=True,
                    help="Index variants")
parser.add_argument("-r", "--reference", action="store", dest="reference", required=True,
                    help="Fasta with reference genome")
parser.add_argument("-c", "--gvcf_input", action="store_true", dest="gvcf_input", default=False,
                    help="Input is gvcf")
parser.add_argument("-g", "--gatk_directory", action="store", dest="gatk_dir", default="",
                    help="Directory with GATK")

args = parser.parse_args()


ValidateVariants4.path = args.gatk_dir
ValidateVariants4.index_vcf(args.reference, args.input_vcf, input_is_gvcf=args.gvcf_input)
