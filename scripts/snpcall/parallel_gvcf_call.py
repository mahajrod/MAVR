#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse

from Tools.GATK import HaplotypeCaller


parser = argparse.ArgumentParser()

parser.add_argument("-b", "--bam_with alignment", action="store", dest="alignment", required=True,
                    help="Bam file with alignment")
parser.add_argument("-o", "--output_dir", action="store", dest="output_dir", required=True,
                    help="Output directory")
parser.add_argument("-p", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-r", "--reference", action="store", dest="reference",
                    help="Fasta with reference genome")
parser.add_argument("-g", "--gatk_directory", action="store", dest="gatk_dir", default="",
                    help="Directory with GATK jar")
parser.add_argument("-t", "--threads", action="store", dest="threads", default=4, type=int,
                    help="Number of threads. Default: 4")
parser.add_argument("-q", "--variant_emit_quality", action="store", dest="emit_quality", default=30, type=int,
                    help="Minimum quality of variant to be emitted")
parser.add_argument("-m", "--memory", action="store", dest="memory", default=10, type=int,
                    help="Maximum memory to use(in gigabytes). Default: 10")

args = parser.parse_args()


HaplotypeCaller.jar_path = args.gatk_dir
#HaplotypeCaller.path = args.path
HaplotypeCaller.max_memory = args.memory
HaplotypeCaller.threads = args.threads

HaplotypeCaller.parallel_gvcf_call(args.reference, args.alignment, args.output_dir, args.output_prefix,
                                   genotyping_mode="DISCOVERY", output_mode="EMIT_VARIANTS_ONLY",
                                   stand_call_conf=args.emit_quality, max_region_length=500000,
                                   max_seqs_per_region=25, length_dict=None, parsing_mode="parse", region_list=None)
