#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse

from Tools.GATK import FastaAlternateReferenceMaker


parser = argparse.ArgumentParser()

parser.add_argument("-v", "--vcf", action="store", dest="vcf", required=True,
                    help="Vcf file")
parser.add_argument("-o", "--output", action="store", dest="output", default="output", required=True,
                    help="New reference")
parser.add_argument("-r", "--reference", action="store", dest="reference", required=True,
                    help="Fasta with reference genome")
parser.add_argument("-g", "--gatk_directory", action="store", dest="gatk_dir", default="", required=True,
                    help="Directory with GATK jar")


args = parser.parse_args()


FastaAlternateReferenceMaker.jar_path = args.gatk_dir
FastaAlternateReferenceMaker.correct_reference(args.reference, args.output, args.vcf)
