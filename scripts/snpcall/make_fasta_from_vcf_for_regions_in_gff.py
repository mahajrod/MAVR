#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.GATK import FastaAlternateReferenceMaker


parser = argparse.ArgumentParser()

parser.add_argument("-v", "--vcf", action="store", dest="vcf", required=True,
                    help="Vcf file")
parser.add_argument("-m", "--mask", action="store", dest="mask",
                    help="Vcf file with masking")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-r", "--reference", action="store", dest="reference", required=True,
                    help="Fasta with reference genome")
parser.add_argument("-f", "--gff", action="store", dest="gff", required=True,
                    help="Gff file with region coordinates")
parser.add_argument("-g", "--gatk_directory", action="store", dest="gatk_dir", default="", required=True,
                    help="Directory with GATK jar")
parser.add_argument("-u", "--use_ambiguous_nucleotides", action="store_true", dest="use_ambiguous_nucleotides",
                    help="Use ambiguous nucleotides for heterozygous positions")

args = parser.parse_args()


FastaAlternateReferenceMaker.jar_path = args.gatk_dir

FastaAlternateReferenceMaker.correct_regions_from_gff(args.reference,
                                                      args.vcf,
                                                      args.gff,
                                                      output_prefix=args.output_prefix,
                                                      feature_type_list=["CDS"],
                                                      unification_key="Parent",
                                                      vcf_with_masking=args.mask,
                                                      override_vcf_by_mask=True,
                                                      use_ambiguous_nuccleotides=args.use_ambiguous_nucleotides)
