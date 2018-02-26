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
parser.add_argument("-s", "--samples_list", action="store", dest="samples_list",
                    help="List of samples to call variants for. Default: all samplews in sample directory")
parser.add_argument("-g", "--gatk_directory", action="store", dest="gatk_dir", default="",
                    help="Directory with GATK jar")

parser.add_argument("-t", "--threads", action="store", dest="threads", default=4, type=int,
                    help="Number of threads. Default: 4")
parser.add_argument("-m", "--memory", action="store", dest="memory", default=10, type=int,
                    help="Maximum memory to use(in gigabytes). Default: 10")
parser.add_argument("-q", "--variant_emit_quality", action="store", dest="emit_quality", default=40, type=int,
                    help="Minimum quality of variant to be emitted")
parser.add_argument("-u", "--variant_call_quality", action="store", dest="call_quality", default=100, type=int,
                    help="Minimum quality of variant to be called. Variant with quality between call and emit will "
                         "be present in vcf file but with filter 'low quality'")
parser.add_argument("-u", "--suffix", action="store", dest="suffix", default="",
                    help="Suffix of the input files. Default: ''")
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


SNPCallPipeline.call_variants(args.sample_dir, args.reference, sample_list=args.sample_list, outdir=args.outdir,
                              suffix=args.suffix, input="alignment",
                              input_filetype="bam", threads=None, mark_duplicates=False,
                              genotyping_mode="DISCOVERY", output_mode="EMIT_VARIANTS_ONLY",
                              stand_emit_conf=args.emit_quality, stand_call_conf=args.call_quality)
