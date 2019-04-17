#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.GATK import UnifiedGenotyper, SelectVariants, VariantFiltration, CombineVariants


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
parser.add_argument("--quality_to_change", action="store", dest="quality_to_change", default=None, type=int,
                    help="Change quality of read alignment with this quality to quality set by "
                         "--quality_to_change_to  on the fly. Default - do not change any read quality")
parser.add_argument("--quality_to_change_to", action="store", dest="quality_to_change_to", default=60, type=int,
                    help="Change quality of read alignment wit quality set by --quality_to_change option "
                         "to this quality on the fly. Works only if --quality_to_change option is set. Default - 60.")


parser.add_argument("--remove_spliced_reads", action="store_true", dest="remove_spliced_reads", default=False,
                    help="Remove spliced reads.")


args = parser.parse_args()


SelectVariants.jar_path = args.gatk_dir
CombineVariants.jar_path = args.gatk_dir
UnifiedGenotyper.jar_path = args.gatk_dir
VariantFiltration.jar_path = args.gatk_dir

UnifiedGenotyper.threads = args.threads

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
                              output_mode="EMIT_VARIANTS_ONLY",
                              discovery_mode="BOTH",
                              output_file=vcf_all,
                              quality_to_change=args.quality_to_change,
                              quality_to_change_to=args.quality_to_change_to,
                              remove_spliced_reads=args.remove_spliced_reads)

SelectVariants.get_indel(args.reference, vcf_all, vcf_indel)
SelectVariants.get_SNP(args.reference, vcf_all, vcf_SNP)

VariantFiltration.filter_bad_SNP(args.reference, vcf_SNP, vcf_filtered_SNP)
VariantFiltration.filter_bad_indel(args.reference, vcf_indel, vcf_filtered_indel)

SelectVariants.remove_entries_with_filters(args.reference, vcf_filtered_SNP, vcf_best_SNP)
SelectVariants.remove_entries_with_filters(args.reference, vcf_filtered_indel, vcf_best_indel)
CombineVariants.combine_from_same_source(args.reference, [vcf_best_SNP, vcf_best_indel], vcf_best_merged)
