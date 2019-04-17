#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Bedtools import Intersect

parser = argparse.ArgumentParser()

parser.add_argument("-a", "--first_vcf", action="store", dest="first_vcf",
                    help="First vcf to compare")
parser.add_argument("-b", "--second_vcf", action="store", dest="second_vcf",
                    help="Second vcf to compare")

parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    help="Prefix of output files")

args = parser.parse_args()

common_vcf = "%s_common.vcf" % args.output_prefix
only_first_vcf = "%s_only_first.vcf" % args.output_prefix
only_second_vcf = "%s_only_second.vcf" % args.output_prefix

Intersect.intersect(args.first_vcf, args.second_vcf, common_vcf, method="-u")
Intersect.intersect(args.first_vcf, args.second_vcf, only_first_vcf, method="-v")
Intersect.intersect(args.second_vcf, args.first_vcf, only_second_vcf, method="-v")