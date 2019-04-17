#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.WGA import LAST


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_fasta_list", action="store", dest="input_fasta_list", required=True,
                    type=LAST.make_list_of_path_to_files_from_string,
                    help="Comma-separated list of input files")
parser.add_argument("-p", "--db_prefix", action="store", dest="db_prefix", required=True,
                    help="Prefix of  LAST database")
parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                    help="Verbose output")
parser.add_argument("-t", "--threads", action="store", dest="threads", default=4, type=int,
                    help="Number of threads. Default: 4")

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
LAST.create_last_db(args.db_prefix,
                    args.input_fasta_list,
                    softmasking=True,
                    seeding_scheme="YASS",
                    verbose=args.verbose,
                    keep_preliminary_masking=True,
                    mask_simple_repeats=True)
