#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.WGA import LAST


parser = argparse.ArgumentParser()

parser.add_argument("-q", "--query_fasta", action="store", dest="query_fasta", required=True,
                    help="Query fasta file")
parser.add_argument("-d", "--db", action="store", dest="db", required=True,
                    help="LAST database")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-f", "--format", action="store", dest="format", default="MAF",
                    help="Format of output file. Allowed: MAF(default), TAB, BlastTab, BlastTab+")
parser.add_argument("-e", "--per_thread_memory", action="store", dest="per_thread_memory",
                    default="4G",
                    help="Per thread memory. Default: 4G")
parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                    help="Verbose output")
parser.add_argument("-t", "--threads", action="store", dest="threads", default=4, type=int,
                    help="Number of threads. Default: 4")
parser.add_argument("-k", "--keep_preliminary_masking", action="store_true",
                    help="Keep preliminary masking. Default: False")
parser.add_argument("-m", "--mask_simple_repeats", action="store_true", dest="mask_simple_repeats",
                    help="Mask simple repeats. Default: False")

parser.add_argument("-g", "--eg2", action="store", dest="eg2", default=None, type=float,
                    help="Maximum EG2 threshold (float). Default: not set")

parser.add_argument("-c", "--cut", action="store", dest="cut", default=None, type=int,
                    help="Lastall -C option. Use it carefully. Default: not set")
"""
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
                    default=4000,  type=int,
                    help="Slurm maximum memmory per cpu in megabytes. Default: 4000")
parser.add_argument("-w", "--slurm_modules_list", action="store", dest="slurm_modules_list", default=[],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of modules to load. Set modules for hmmer and python")
"""

args = parser.parse_args()

LAST.threads = args.threads
LAST.lastal(args.db, args.query_fasta, args.output_prefix, verbose=args.verbose,
            keep_preliminary_masking=args.keep_preliminary_masking,
            mask_simple_repeats=args.mask_simple_repeats,
            output_format=args.format, per_thread_memory=args.per_thread_memory,
            eg2_threshold=args.eg2, discard_limit=args.cut)
