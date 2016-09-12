#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse

from Routines import FileRoutines
from Tools.GATK import VariantFiltration

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_vcf", action="store", dest="input_vcf", required=True,
                    help="Input vcf file")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-r", "--reference", action="store", dest="reference", required=True,
                    help="Fasta with reference genome")
parser.add_argument("-g", "--gatk_directory", action="store", dest="gatk_dir", default="",
                    help="Directory with GATK jar")
parser.add_argument("--snp_filter_name", action="store", dest="snp_filter_name", type=str,
                    default="ambiguous_snp", help="SNP filter name")
parser.add_argument("--snp_QD", action="store", dest="snp_QD", type=float, default=2.0,
                    help="SNP QD threshold. Default -  2.0")
parser.add_argument("--snp_FS", action="store", dest="snp_FS", type=float, default=60.0,
                    help="SNP FS threshold. Default -   60.0")
parser.add_argument("--snp_MQ", action="store", dest="snp_MQ", type=float, default=40.0,
                    help="SNP MQ threshold. Default -  40.0")
parser.add_argument("--snp_HaplotypeScore", action="store", dest="snp_HaplotypeScore", type=float, default=13.0,
                    help="SNP HaplotypeScore threshold. Default -  13.0")
parser.add_argument("--snp_MappingQualityRankSum", action="store", dest="snp_MappingQualityRankSum", type=float,
                    default=-12.5, help="SNP MappingQualityRankSum threshold. Default -")
parser.add_argument("--snp_ReadPosRankSum", action="store", dest="snp_ReadPosRankSum", type=float, default=-8.0,
                    help="SNP ReadPosRankSum threshold. Default -   -8.0")
parser.add_argument("--indel_filter_name", action="store", dest="indel_filter_name", type=str,
                    default="ambiguous_indel", help="Indel filter name")
parser.add_argument("--indel_QD", action="store", dest="indel_QD", type=float, default=2.0,
                    help="Indel QD threshold. Default - 2.0")
parser.add_argument("--indel_ReadPosRankSum", action="store", dest="indel_ReadPosRankSum", type=float, default=-20.0,
                    help="Indel ReadPosRankSum threshold. Default -   -20.0")
parser.add_argument("--indel_InbreedingCoeff", action="store", dest="indel_InbreedingCoeff", type=float, default=-0.8,
                    help="Indel InbreedingCoeff threshold. Default -   -0.8")
parser.add_argument("--indel_FS", action="store", dest="indel_FS", type=float, default=200.0,
                    help="Indel FS threshold. Default - 200.0")

args = parser.parse_args()

VariantFiltration.jar_path = FileRoutines.check_path(args.gatk_dir)

VariantFiltration.filter_bad_variants(args.reference, args.input_vcf, args.output_prefix,
                                      snp_filter_name=args.snp_filter_name, snp_QD=args.snp_QD,
                                      snp_FS=args.snp_FS, snp_MQ=args.snp_MQ,
                                      snp_HaplotypeScore=args.snp_HaplotypeScore,
                                      snp_MappingQualityRankSum=-args.snp_MappingQualityRankSum,
                                      snp_ReadPosRankSum=args.snp_ReadPosRankSum,
                                      indel_filter_name=args.indel_filter_name, indel_QD=args.indel_QD,
                                      indel_ReadPosRankSum=args.indel_ReadPosRankSum,
                                      indel_InbreedingCoeff=args.indel_InbreedingCoeff, indel_FS=args.indel_FS)

