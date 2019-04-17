#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse
from RouToolPa.Tools.Bedtools import Intersect



parser = argparse.ArgumentParser()

parser.add_argument("-a", "--first_vcf", action="store", dest="first_vcf", required=True,
                    help="First vcf file")
parser.add_argument("-b", "--second_vcf", action="store", dest="second_vcf", required=True,
                    help="Second vcf file")
parser.add_argument("-o", "--output_vcf", action="store", dest="output_vcf", required=True,
                    help="Output vcf file")

args = parser.parse_args()


with open(args.first_vcf, "r") as in_fd:
    with open("header.tmp", "w") as h_fd:
        for line in in_fd:
            if line[0] != "#":
                break
            h_fd.write(line)

Intersect.intersect(args.first_vcf, args.second_vcf, "data.vcf", method="-wa")  # entries from a which intersect with entries from b
os.system("cat %s %s > %s" % ("header.tmp", "data.vcf", args.output_vcf))

os.remove("header.tmp")
os.remove("data.vcf")