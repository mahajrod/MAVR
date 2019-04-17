#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse
from Pipelines.TenX import TenXAlignmentPipeline
from RouToolPa.Routines import FileRoutines



parser = argparse.ArgumentParser()

parser.add_argument("-d", "--sample_directory", action="store", dest="samples_dir", required=True,
                    type=lambda s: FileRoutines.check_path(os.path.abspath(s)),
                    help="Directory with samples")
parser.add_argument("-s", "--samples", action="store", dest="samples", type=lambda s: s.split(","),
                    help="Comma-separated list of subdirectories(one per sample) to handle. "
                         "If not set all subdirectories will be considered as containing samples")
parser.add_argument("-o", "--output_dir", action="store", dest="output_dir",
                    type=lambda s: FileRoutines.check_path(os.path.abspath(s)),
                    default="./", help="Directory to write output. Default: current directory")
parser.add_argument("-t", "--threads", action="store", dest="threads", default=1, type=int,
                    help="Number of threads to use in Trimmomatic. Default - 1.")
parser.add_argument("-r", "--reference", action="store", dest="reference",
                    required=True, help="File with genome reference")
parser.add_argument("-v", "--variant_calling_mode", action="store", dest="variant_calling_mode",
                    help="Variant calling mode. Allowed: freebayes(default), gatk(-g/--gatk_jar option is required),"
                         " disable (phase only)")
parser.add_argument("-g", "--gatk_jar", action="store", dest="gatk_jar",
                    help="Path to GATK jar. Required only if variant calling with GATK was chosen.")
parser.add_argument("-m", "--max_memory", action="store", dest="max_memory", type=int, default=10,
                    help="Maximum memory usage allowed (in gigabytes). Default - 10.")
parser.add_argument("-p", "--precalled_vcf", action="store", dest="precalled_vcf",
                    help="Use this precalled vcf for phasing")
parser.add_argument("-x", "--sex", action="store", dest="sex",
                    help="Sex of the sample. Allowed: m, f, male or female. Auto detection if not set")
parser.add_argument("-l", "--longranger_dir", action="store", dest="longranger_dir",
                    help="Path to directory with LongRanger binary")
parser.add_argument("-a", "--use_somatic_sv_caller", action="store_true", dest="use_somatic_sv_caller", default=False,
                    help="Use somatic structural variant caller")
parser.add_argument("-c", "--variant_calling_only", action="store_true", dest="variant_calling_only", default=False,
                    help="Call variants only. No structural variant calling of phasing")

"""
parser.add_argument("-x", "--general_stat_file", action="store", dest="general_stat_file", required=True,
                    help="File to write general statistics about filtration")
"""

args = parser.parse_args()

if args.gatk_jar:
    args.variant_calling_mode = "gatk"


TenXAlignmentPipeline.threads = args.threads
TenXAlignmentPipeline.max_memory = args.max_memory
TenXAlignmentPipeline.longranger_dir = args.longranger_dir
TenXAlignmentPipeline.align_and_call(args.samples_dir, args.output_dir, args.reference, samples_to_handle=args.samples,
                                     variant_calling_mode=args.variant_calling_mode, gatk_jar_path=args.gatk_jar,
                                     use_somatic_sv_caller=args.use_somatic_sv_caller, precalled_vcf=args.precalled_vcf,
                                     sample_sex=args.sex, variant_calling_only=args.variant_calling_only)

