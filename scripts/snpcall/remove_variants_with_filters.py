#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.GATK import SelectVariants
from RouToolPa.Routines import FileRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_vcf", action="store", dest="input_vcf", required=True,
                    help="Input vcf file")
parser.add_argument("-o", "--output_vcf", action="store", dest="output_vcf", required=True,
                    help="Output vcf file")
parser.add_argument("-r", "--reference", action="store", dest="reference", required=True,
                    help="Fasta with reference genome")
parser.add_argument("-g", "--gatk_directory", action="store", dest="gatk_dir", default="",
                    help="Directory with GATK jar")
args = parser.parse_args()

SelectVariants.jar_path = FileRoutines.check_path(args.gatk_dir)
SelectVariants.remove_entries_with_filters(args.reference, args.input_vcf, args.output_vcf)
