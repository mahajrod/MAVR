#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse

from Pipelines import SNPCallPipeline


parser = argparse.ArgumentParser()

parser.add_argument("-d", "--sample_dir", action="store", dest="sample_dir", required=True,
                    help="Directory with samples data")
parser.add_argument("-r", "--reference", action="store", dest="reference", required=True,
                    help="File with reference genome")
parser.add_argument("-o", "--outdir", action="store", dest="outdir", default="./",
                    help="Output directory. Default: current directory")
parser.add_argument("-e", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of merged files")
parser.add_argument("-s", "--sample_list", action="store", dest="sample_list",
                    help="List of samples to call variants for. Default: all samplews in sample directory")
parser.add_argument("-g", "--gatk_directory", action="store", dest="gatk_dir", default="",
                    help="Directory with GATK jar")
parser.add_argument("-p", "--picard_directory", action="store", dest="picard_dir", default="",
                    help="Directory with PICARD jar")

parser.add_argument("-t", "--threads", action="store", dest="threads", default=4, type=int,
                    help="Number of threads. Default: 4")
parser.add_argument("-m", "--memory", action="store", dest="memory", default=10, type=int,
                    help="Maximum memory to use(in gigabytes). Default: 10")
parser.add_argument("-u", "--variant_call_quality", action="store", dest="call_quality", default=30, type=int,
                    help="Minimum quality of variant to be called.") # Variant with quality between call and emit will ""be present in vcf file but with filter 'low quality'")
parser.add_argument("-x", "--suffix", action="store", dest="suffix", default="",
                    help="Suffix of the input files. Default: ''")

parser.add_argument("-y", "--regions_to_include", action="store", dest="regions_to_include",
                    help="File with regions to include in analysis")
parser.add_argument("-z", "--regions_to_exclude", action="store", dest="regions_to_exclude",
                    help="File with regions to exclude in analysis")

"""
parser.add_argument("--quality_to_change", action="store", dest="quality_to_change", default=None, type=int,
                    help="Change quality of read alignment with this quality to quality set by "
                         "--quality_to_change_to  on the fly. Default - do not change any read quality")
parser.add_argument("--quality_to_change_to", action="store", dest="quality_to_change_to", default=60, type=int,
                    help="Change quality of read alignment wit quality set by --quality_to_change option "
                         "to this quality on the fly. Works only if --quality_to_change option is set. Default - 60.")


parser.add_argument("--remove_spliced_reads", action="store_true", dest="remove_spliced_reads", default=False,
                    help="Remove spliced reads.")
"""

args = parser.parse_args()

SNPCallPipeline.threads = args.threads
SNPCallPipeline.max_memory = args.memory
SNPCallPipeline.GATK_dir = args.gatk_dir
SNPCallPipeline.Picard_dir = args.picard_dir

SNPCallPipeline.call_variants(args.sample_dir, args.reference, args.output_prefix, sample_list=args.sample_list, outdir=args.outdir,
                              suffix=args.suffix, input="alignment",
                              input_filetype="bam", threads=None, mark_duplicates=False,
                              genotyping_mode="DISCOVERY", output_mode="EMIT_VARIANTS_ONLY",
                              stand_call_conf=args.call_quality,
                              skip_base_score_recalibration=False,
                              iteration_number=3, SNP_QD=2.0, SNP_FS=30.0, SNP_MQ=40.0, SNP_MappingQualityRankSum=-12.5,
                              SNP_ReadPosRankSum=-8.0, indel_QD=2.0, indel_ReadPosRankSum=-20.0, indel_FS=200.0,
                              SNP_filter_name="ambiguous_snp", indel_filter_name="ambiguous_indel",
                              analyze_covariates=True,
                              include_region_id_file=args.regions_to_include,
                              exclude_region_id_file=args.regions_to_exclude)
