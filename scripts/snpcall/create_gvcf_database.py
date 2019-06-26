#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.GATK4 import GenomicsDBImport4

parser = argparse.ArgumentParser()

parser.add_argument("-o", "--output_dbi_dir", action="store", dest="output_dbi_dir", required=True,
                    help="Output directory to write database")
parser.add_argument("-g", "--gatk_directory", action="store", dest="gatk_dir", default="",
                    help="Directory with GATK")
parser.add_argument("-i", "--gvcf_list", action="store", dest="gvcf_list",
                    type=GenomicsDBImport4.make_list_of_path_to_files_from_string,
                    help="Comma-separated list of gvcf files to include in database",
                    required=True)
parser.add_argument("-l", "--interval_list", action="store", dest="interval_list", required=True,
                    help="Comma-separated list of intervals to import.")
parser.add_argument("-s", "--extension_list", action="store", dest="extension_list", default=["g.vcf",],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of extension of GVCF files. Default: g.vcf")
parser.add_argument("-m", "--memory", action="store", dest="memory", default="10g",
                    help="Maximum memory to use. Default: 10g")

"""


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
"""
args = parser.parse_args()

GenomicsDBImport4.path = args.gatk_dir
GenomicsDBImport4.threads = 1
GenomicsDBImport4.max_memory = args.memory

GenomicsDBImport4.create_db(args.gvcf_list,
                            args.output_dbi_dir,
                            interval_list=args.interval_list,
                            extension_list=args.extension_list)
