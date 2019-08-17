#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.GATK import VariantFiltration
from RouToolPa.Tools.GATK4 import VariantFiltration4
from RouToolPa.Routines import FileRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_vcf", action="store", dest="input_vcf", required=True,
                    help="Input vcf file")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-r", "--reference", action="store", dest="reference", required=True,
                    help="Fasta with reference genome")
parser.add_argument("-c", "--combine", action="store_true", dest="combine",
                    help="Combine indel and SNP vcfs after filtration")
parser.add_argument("-s", "--sequence_dict_file", action="store", dest="sequence_dict_file",
                    help="Sequence dictionary file. Required for sorting of merged files according to the dict")
parser.add_argument("-g", "--gatk_directory", action="store", dest="gatk_dir", default="",
                    help="Directory with GATK jar")
parser.add_argument("-p", "--picard_directory", action="store", dest="picard_dir", default="",
                    help="Directory with Picard jar")

parser.add_argument("-m", "--memory", action="store", dest="memory", default="10000", type=lambda s: s + "m",
                    help="Maximum memory to use in megabytes. Default: 10000")

parser.add_argument("--snp_filter_name", action="store", dest="snp_filter_name", type=str,
                    default="ambiguous_snp", help="SNP filter name")
parser.add_argument("--snp_QD", action="store", dest="snp_QD", type=float, default=2.0,
                    help="SNP QD threshold. Default -  2.0")
parser.add_argument("--snp_FS", action="store", dest="snp_FS", type=float, default=60.0,
                    help="SNP FS threshold. Default -   60.0")
parser.add_argument("--snp_MQ", action="store", dest="snp_MQ", type=float, default=40.0,
                    help="SNP MQ threshold. Default -  40.0")
#parser.add_argument("--snp_HaplotypeScore", action="store", dest="snp_HaplotypeScore", type=float, default=13.0,
#                    help="SNP HaplotypeScore threshold. Default -  13.0")
parser.add_argument("--snp_MappingQualityRankSum", action="store", dest="snp_MappingQualityRankSum", type=float,
                    default=-12.5, help="SNP MappingQualityRankSum threshold. Default - -12.5")
parser.add_argument("--snp_ReadPosRankSum", action="store", dest="snp_ReadPosRankSum", type=float, default=-8.0,
                    help="SNP ReadPosRankSum threshold. Default -   -8.0")
parser.add_argument("--indel_filter_name", action="store", dest="indel_filter_name", type=str,
                    default="ambiguous_indel", help="Indel filter name")
parser.add_argument("--indel_QD", action="store", dest="indel_QD", type=float, default=2.0,
                    help="Indel QD threshold. Default - 2.0")
parser.add_argument("--indel_ReadPosRankSum", action="store", dest="indel_ReadPosRankSum", type=float, default=-20.0,
                    help="Indel ReadPosRankSum threshold. Default -   -20.0")
#parser.add_argument("--indel_InbreedingCoeff", action="store", dest="indel_InbreedingCoeff", type=float, default=-0.8,
#                    help="Indel InbreedingCoeff threshold. Default -   -0.8")
parser.add_argument("--indel_FS", action="store", dest="indel_FS", type=float, default=200.0,
                    help="Indel FS threshold. Default - 200.0")
parser.add_argument("-v", "--gatk_version", action="store", dest="gatk_version", default="4",
                    help="Major version of GATK. Allowed: 4(default), 3 ")

args = parser.parse_args()

filtration_tool = VariantFiltration4 if args.gatk_version == "4" else VariantFiltration

filtration_tool.path = FileRoutines.check_path(args.gatk_dir)

filtration_tool.filter_bad_variants(args.reference, args.input_vcf, args.output_prefix,
                                      snp_filter_name=args.snp_filter_name, snp_QD=args.snp_QD,
                                      snp_FS=args.snp_FS, snp_MQ=args.snp_MQ,
                                      #snp_HaplotypeScore=args.snp_HaplotypeScore,
                                      snp_MappingQualityRankSum=args.snp_MappingQualityRankSum,
                                      snp_ReadPosRankSum=args.snp_ReadPosRankSum,
                                      indel_filter_name=args.indel_filter_name, indel_QD=args.indel_QD,
                                      indel_ReadPosRankSum=args.indel_ReadPosRankSum,
                                      indel_FS=args.indel_FS,
                                      combine_vcf=args.combine, sequence_dict_file=args.sequence_dict_file,
                                      picard_memory=args.memory, picard_dir=args.picard_dir)
                                      #indel_InbreedingCoeff=args.indel_InbreedingCoeff, )

