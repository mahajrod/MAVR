#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.GATK4 import SortVcf4

parser = argparse.ArgumentParser()

parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output sorted file")
parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Comma-separated list of input files")

parser.add_argument("-s", "--sequence_dict", action="store", dest="sequence_dict",
                    help="File with sequence dict corresponding to the reference genome")
parser.add_argument("-g", "--gatk_directory", action="store", dest="gatk_dir", default="",
                    help="Directory with GATK")

parser.add_argument("-m", "--memory", action="store", dest="memory", default="10g",
                    help="Maximum memory to use. Default: 10g")

parser.add_argument("-d", "--handling_mode", action="store", dest="handling_mode", default="local",
                    help="Handling mode. Allowed: local(default), slurm")
parser.add_argument("-j", "--slurm_job_name", action="store", dest="slurm_job_name", default="JOB",
                    help="Slurm job name. Default: JOB")

parser.add_argument("-y", "--slurm_log_prefix", action="store", dest="slurm_log_prefix",
                    help="Slurm log prefix. ")
parser.add_argument("-l", "--slurm_error_log_prefix", action="store", dest="slurm_error_log_prefix",
                    help="Slurm error log prefix")

parser.add_argument("-a", "--slurm_max_running_time", action="store", dest="slurm_max_running_time", default="100:00:00",
                    help="Slurm max running time in hh:mm:ss format. Default: 100:00:00")

parser.add_argument("-w", "--slurm_modules_list", action="store", dest="slurm_modules_list", default=[],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of modules to load. Set modules for hmmer and python")
parser.add_argument("-e", "--tmp_dir", action="store", dest="tmp_dir",
                    help="Temporary directory")

args = parser.parse_args()

SortVcf4.path = args.gatk_dir
SortVcf4.threads = 1
SortVcf4.max_memory = args.memory
SortVcf4.tmp_dir = args.tmp_dir

SortVcf4.sort_vcf(args.input, args.output, seq_dict=args.sequence_dict,
                  handling_mode=args.handling_mode,
                  max_memory_per_node=args.memory,
                  max_running_time=args.slurm_max_running_time,
                  job_name=args.slurm_job_name,
                  log_prefix=args.slurm_log_prefix,
                  error_log_prefix=args.slurm_error_log_prefix,
                  modules_list=args.slurm_modules_list,
                  environment_variables_dict=None)


