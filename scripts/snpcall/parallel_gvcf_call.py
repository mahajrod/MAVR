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
                    help="Minimum quality of variant to be emitted. Default: 30")
parser.add_argument("-m", "--memory", action="store", dest="memory", default=10, type=int,
                    help="Maximum memory to use(in gigabytes). Default: 10")
parser.add_argument("-x", "--max_seq_per_region", action="store", dest="max_seq_per_region", default=10, type=int,
                    help="Maximum sequences per region. Default: 10")
parser.add_argument("-l", "--max_region_len", action="store", dest="max_region_len", default=50000, type=int,
                    help="Maximum region length. Default: 50000")

parser.add_argument("-d", "--handling_mode", action="store", dest="handling_mode", default="local",
                    help="Handling mode. Allowed: local(default), slurm")
parser.add_argument("-j", "--slurm_job_name", action="store", dest="slurm_job_name", default="JOB",
                    help="Slurm job name. Default: JOB")

parser.add_argument("-y", "--slurm_log_prefix", action="store", dest="slurm_log_prefix",
                    help="Slurm log prefix. ")
parser.add_argument("-e", "--slurm_error_log_prefix", action="store", dest="slurm_error_log_prefix",
                    help="Slurm error log prefix")
parser.add_argument("-z", "--slurm_max_running_jobs", action="store", dest="slurm_max_running_jobs",
                    default=300, type=int,
                    help="Slurm max running jobs. Default: 300")
parser.add_argument("-a", "--slurm_max_running_time", action="store", dest="slurm_max_running_time", default="100:00:00",
                    help="Slurm max running time in hh:mm:ss format. Default: 100:00:00")

parser.add_argument("-u", "--slurm_max_memmory_per_cpu", action="store", dest="slurm_max_memmory_per_cpu",
                    default=4000, type=int,
                    help="Slurm maximum memmory per cpu in megabytes. Default: 4000")
parser.add_argument("-q", "--slurm_modules_list", action="store", dest="slurm_modules_list", default=[],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of modules to load. Set modules for hmmer and python")

args = parser.parse_args()


HaplotypeCaller.jar_path = args.gatk_dir
#HaplotypeCaller.path = args.path
HaplotypeCaller.max_memory = args.memory
HaplotypeCaller.threads = args.threads

HaplotypeCaller.parallel_gvcf_call(args.reference, args.alignment, args.output_dir, args.output_prefix,
                                   "%s.combined.g.vcf" % args.output_prefix,
                                   genotyping_mode="DISCOVERY",
                                   stand_call_conf=args.emit_quality, max_region_length=args.max_region_len,
                                   max_seqs_per_region=args.max_seq_per_region, length_dict=None,
                                   parsing_mode="parse", region_list=None,
                                   handling_mode=args.handling_mode,
                                   job_name=args.slurm_job_name,
                                   log_prefix=args.slurm_log_prefix,
                                   error_log_prefix=args.slurm_error_log_prefix,
                                   max_running_jobs=args.slurm_max_running_jobs,
                                   max_running_time=args.slurm_max_running_time,
                                   max_memmory_per_cpu=args.slurm_max_memmory_per_cpu,)

"""
parallel_gvcf_call(self, reference, alignment, output_dir, output_prefix, output,
                           genotyping_mode="DISCOVERY",
                           stand_call_conf=30, max_region_length=1000000, max_seqs_per_region=100,
                           length_dict=None, parsing_mode="parse", region_list=None, remove_intermediate_files=False,
                           gvcf_extension_list=["g.vcf", ],)
                           """