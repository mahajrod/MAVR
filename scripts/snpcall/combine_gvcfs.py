#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.GATK4 import CombineGVCFs4

parser = argparse.ArgumentParser()

parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output directory to write database")

parser.add_argument("-i", "--gvcf_list", action="store", dest="gvcf_list",
                    type=CombineGVCFs4.make_list_of_path_to_files_from_string,
                    help="Comma-separated list of gvcf files to include in database",
                    required=True)
parser.add_argument("-r", "--reference", action="store", dest="reference", required=True,
                    help="Fasta with reference genome")
parser.add_argument("-s", "--extension_list", action="store", dest="extension_list", default=["g.vcf",],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of extension of GVCF files. Default: g.vcf")
parser.add_argument("-m", "--memory", action="store", dest="memory", default="10g",
                    help="Maximum memory to use. Default: 10g")
parser.add_argument("-g", "--gatk_directory", action="store", dest="gatk_dir", default="",
                    help="Directory with GATK")
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

CombineGVCFs4.path = args.gatk_dir
CombineGVCFs4.threads = 1
CombineGVCFs4.max_memory = args.memory

CombineGVCFs4.combine(args.reference,
                      args.gvcf_list,
                      args.output,
                      extension_list=args.extension_list)
