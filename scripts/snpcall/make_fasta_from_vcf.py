#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.GATK import FastaAlternateReferenceMaker


parser = argparse.ArgumentParser()

parser.add_argument("-v", "--vcf", action="store", dest="vcf", required=True,
                    help="Vcf file")
parser.add_argument("-o", "--output", action="store", dest="output", default="output", required=True,
                    help="New reference")
parser.add_argument("-r", "--reference", action="store", dest="reference", required=True,
                    help="Fasta with reference genome")
parser.add_argument("-g", "--gatk_directory", action="store", dest="gatk_dir", default="", required=True,
                    help="Directory with GATK jar")
parser.add_argument("-m", "--memory", action="store", dest="memory", default="10g",
                    help="Maximum memory to use. Default: 10g")
parser.add_argument("-a", "--mask", action="store", dest="mask",
                    help="Vcf file with masking")
parser.add_argument("-u", "--use_ambiguous_nucleotides", action="store_true", dest="use_ambiguous_nucleotides",
                    help="Use ambiguous nucleotides for heterozygous positions")
args = parser.parse_args()


FastaAlternateReferenceMaker.jar_path = args.gatk_dir
FastaAlternateReferenceMaker.max_memory = args.memory

FastaAlternateReferenceMaker.correct_reference(args.reference, args.output, args.vcf,
                                               vcf_with_masking=args.mask,
                                               use_ambiguous_nuccleotides=args.use_ambiguous_nucleotides)
