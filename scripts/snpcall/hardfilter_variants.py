#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse

from Pipelines import SNPCallPipeline


parser = argparse.ArgumentParser()

parser.add_argument("-r", "--reference", action="store", dest="reference", required=True,
                    help="File with reference genome")
parser.add_argument("-i", "--input_vcf", action="store", dest="input_vcf", required=True,
                    help="Input vcf file")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

parser.add_argument("--snp_filter_name", action="store", dest="snp_filter_name", type=str,
                    default="ambiguous_snp", help="SNP filter name. Default: ambiguous_snp")
parser.add_argument("--snp_QD", action="store", dest="snp_QD", type=float, default=2.0,
                    help="SNP QD threshold. Default -  2.0")
parser.add_argument("--snp_FS", action="store", dest="snp_FS", type=float, default=30.0,
                    help="SNP FS threshold. Default -   30.0")
parser.add_argument("--snp_MQ", action="store", dest="snp_MQ", type=float, default=40.0,
                    help="SNP MQ threshold. Default -  40.0")
parser.add_argument("--snp_MappingQualityRankSum", action="store", dest="snp_MappingQualityRankSum", type=float,
                    default=-12.5, help="SNP MappingQualityRankSum threshold. Default - -12.5")
parser.add_argument("--snp_ReadPosRankSum", action="store", dest="snp_ReadPosRankSum", type=float, default=-8.0,
                    help="SNP ReadPosRankSum threshold. Default -   -8.0")

parser.add_argument("--indel_filter_name", action="store", dest="indel_filter_name", type=str,
                    default="ambiguous_indel", help="Indel filter name. Default: ambiguous indel")
parser.add_argument("--indel_QD", action="store", dest="indel_QD", type=float, default=2.0,
                    help="Indel QD threshold. Default - 2.0")
parser.add_argument("--indel_ReadPosRankSum", action="store", dest="indel_ReadPosRankSum", type=float, default=-20.0,
                    help="Indel ReadPosRankSum threshold. Default -   -20.0")
parser.add_argument("--indel_FS", action="store", dest="indel_FS", type=float, default=200.0,
                    help="Indel FS threshold. Default - 200.0")


parser.add_argument("-g", "--gatk_directory", action="store", dest="gatk_dir", default="",
                    help="Directory with GATK jar")
parser.add_argument("-t", "--threads", action="store", dest="threads", default=4, type=int,
                    help="Number of threads. Default: 4")
parser.add_argument("-m", "--memory", action="store", dest="memory", default=10, type=int,
                    help="Maximum memory to use(in gigabytes). Default: 10")

"""
parser.add_argument("-y", "--regions_to_include", action="store", dest="regions_to_include",
                    help="File with regions to include in analysis")
parser.add_argument("-z", "--regions_to_exclude", action="store", dest="regions_to_exclude",
                    help="File with regions to exclude in analysis")
"""


args = parser.parse_args()

SNPCallPipeline.threads = args.threads
SNPCallPipeline.max_memory = args.memory
SNPCallPipeline.GATK_dir = args.gatk_dir


SNPCallPipeline.hardfilter_variants(args.reference, args.input_vcf, args.output_prefix, SNP_QD=args.snp_QD,
                                    SNP_FS=args.snp_FS, SNP_MQ=args.snp_MQ,
                                    SNP_MappingQualityRankSum=args.snp_MappingQualityRankSum,
                                    SNP_ReadPosRankSum=args.snp_ReadPosRankSum,
                                    indel_QD=args.indel_QD, indel_ReadPosRankSum=args.indel_ReadPosRankSum,
                                    indel_FS=args.indel_FS,
                                    SNP_filter_name=args.snp_filter_name,
                                    indel_filter_name=args.indel_filter_name, threads=None)
