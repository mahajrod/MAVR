#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse

from Tools.GATK import UnifiedGenotyper, SelectVariants, VariantFiltration,CombineVariants


parser = argparse.ArgumentParser()

parser.add_argument("-b", "--bam_with alignment", action="store", dest="alignment",
                    help="Bam file with alignment")
parser.add_argument("-p", "--prefix", action="store", dest="prefix", default="alignment",
                    help="Prefix of output files")
parser.add_argument("-r", "--reference", action="store", dest="reference",
                    help="Fasta with reference genome")
parser.add_argument("-g", "--gatk_directory", action="store", dest="gatk_dir", default="",
                    help="Directory with GATK jar")
parser.add_argument("-t", "--threads", action="store", dest="threads", default=4, type=int,
                    help="Number of threads. Default: 4")
parser.add_argument("-q", "--variant_emit_quality", action="store", dest="emit_quality", default=40, type=int,
                    help="Minimum quality of variant to be emitted")
parser.add_argument("-u", "--variant_call_quality", action="store", dest="call_quality", default=100, type=int,
                    help="Minimum quality of variant to be called. Variant with quality between call and emit will "
                         "be present in vcf file but with filter 'low quality'")

args = parser.parse_args()

rmdup_sorted_filtered_alignment = "%s_final.bam" % args.prefix

vcf_all = "%s_GATK_raw.vcf" % args.prefix
vcf_indel = "%s_GATK_raw_indel.vcf" % args.prefix
vcf_SNP = "%s_GATK_raw_SNP.vcf" % args.prefix
vcf_filtered_indel = "%s_GATK_filtered_indel.vcf" % args.prefix
vcf_filtered_SNP = "%s_GATK_filtered_SNP.vcf" % args.prefix
vcf_best_indel = "%s_GATK_best_indel.vcf" % args.prefix
vcf_best_SNP = "%s_GATK_best_SNP.vcf" % args.prefix
vcf_best_merged = "%s_GATK_best_merged.vcf" % args.prefix

UnifiedGenotyper.variant_call(args.alignment,
                              args.reference,
                              stand_emit_conf=args.emit_quality,
                              stand_call_conf=args.call_quality,
                              GATK_dir=args.gatk_dir,
                              num_of_threads=args.threads,
                              output_mode="EMIT_VARIANTS_ONLY",
                              discovery_mode="BOTH",
                              output_file=vcf_all)

SelectVariants.get_indel(args.gatk_dir, args.reference, vcf_all, vcf_indel)
SelectVariants.get_SNP(args.gatk_dir, args.reference, vcf_all, vcf_SNP)

VariantFiltration.filter_bad_SNP(args.gatk_dir, args.reference, vcf_SNP, vcf_filtered_SNP)
VariantFiltration.filter_bad_indel(args.gatk_dir, args.reference, vcf_indel, vcf_filtered_indel)
SelectVariants.remove_filtered(args.gatk_dir, args.reference, vcf_filtered_SNP, vcf_best_SNP)
SelectVariants.remove_filtered(args.gatk_dir, args.reference, vcf_filtered_indel, vcf_best_indel)
CombineVariants.combine_from_same_source(args.gatk_dir, args.reference, [vcf_best_SNP, vcf_best_indel], vcf_best_merged)