#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse

from Tools.GATK import GenotypeGVCFs

parser = argparse.ArgumentParser()

parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

parser.add_argument("-r", "--reference", action="store", dest="reference", required=True,
                    help="Fasta with reference genome")
parser.add_argument("-g", "--gatk_directory", action="store", dest="gatk_dir", default="",
                    help="Directory with GATK jar")
parser.add_argument("-i", "--gvcf_list", action="store", dest="gvcf_list",
                    type=GenotypeGVCFs.make_list_of_path_to_files_from_string,
                    help="Comma-separated list of gvcf files to genotype",  required=True,)

parser.add_argument("-e", "--extension_list", action="store", dest="extension_list", default=["g.vcf",],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of extension of GVCF files. Default: g.vcf")

parser.add_argument("-x", "--max_alternate_alleles", action="store", dest="max_alternate_alleles", type=int,
                    help="Maximum number of alternative allels. Default: GATK default")

parser.add_argument("-m", "--memory", action="store", dest="memory", default="10000", type=lambda s: s + "m",
                    help="Maximum memory to use in megabytes. Default: 10000")

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
args = parser.parse_args()

GenotypeGVCFs.jar_path = args.gatk_dir
GenotypeGVCFs.threads = 1
GenotypeGVCFs.max_memory = args.memory

GenotypeGVCFs.genotype(args.reference,
                       args.gvcf_list,
                       args.output_prefix,
                       extension_list=args.extension_list,
                       handling_mode=args.handling_mode,
                       max_memory_per_node=args.memory,
                       max_running_time=args.slurm_max_running_time,
                       job_name=args.slurm_job_name,
                       log_prefix=args.slurm_log_prefix,
                       error_log_prefix=args.slurm_error_log_prefix,
                       modules_list=args.slurm_modules_list,
                       environment_variables_dict=None,
                       max_alternate_alleles=args.max_alternate_alleles)


