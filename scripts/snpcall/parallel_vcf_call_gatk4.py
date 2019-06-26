#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.GATK import HaplotypeCaller
from RouToolPa.Tools.GATK4 import HaplotypeCaller4

parser = argparse.ArgumentParser()

parser.add_argument("-b", "--alignment_list", action="store", dest="alignment", required=True,
                    type=HaplotypeCaller4.make_list_of_path_to_files_from_string,
                    help="Comma-separated list of bams")
parser.add_argument("-o", "--output_dir", action="store", dest="output_dir", required=True,
                    help="Output directory")
parser.add_argument("-p", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-r", "--reference", action="store", dest="reference", required=True,
                    help="Fasta with reference genome")
parser.add_argument("-v", "--gatk_version", action="store", dest="gatk_version", default="4",
                    help="Major version of GATK. Allowed: 4(default), 3 ")
parser.add_argument("-g", "--gatk_directory", action="store", dest="gatk_dir", default="",
                    help="Directory with GATK")
parser.add_argument("-t", "--threads", action="store", dest="threads", default=4, type=int,
                    help="Number of threads. Default: 4")
parser.add_argument("-q", "--variant_emit_quality", action="store", dest="emit_quality", default=30, type=int,
                    help="Minimum quality of variant to be emitted. Default: 30")
parser.add_argument("-m", "--memory", action="store", dest="memory", default="10g",
                    help="Maximum memory to use. Default: 10g")
parser.add_argument("-x", "--max_seq_per_region", action="store", dest="max_seq_per_region", default=100, type=int,
                    help="Maximum sequences per region. Default: 100")
parser.add_argument("-l", "--max_region_len", action="store", dest="max_region_len", default=2000000, type=int,
                    help="Maximum region length. Default: 2000000")
parser.add_argument("-k", "--black_list_scaffold_id_file", action="store", dest="black_list_scaffold_id_file",
                    help="Id file with scaffolds from black list")

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
parser.add_argument("-a", "--slurm_max_running_time", action="store", dest="slurm_max_running_time",
                    default="100:00:00",
                    help="Slurm max running time in hh:mm:ss format. Default: 100:00:00")

parser.add_argument("-u", "--slurm_max_memmory_per_cpu", action="store", dest="slurm_max_memmory_per_cpu",
                    default=4000,  type=int,
                    help="Slurm maximum memmory per cpu in megabytes. Default: 4000")
parser.add_argument("-w", "--slurm_modules_list", action="store", dest="slurm_modules_list", default=[],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of modules to load. Set modules for hmmer and python")

args = parser.parse_args()

HaplotypeCaller4.max_memory = args.memory
HaplotypeCaller4.threads = args.threads
HaplotypeCaller4.path = args.gatk_dir
HaplotypeCaller4.parallel_call(args.reference, args.alignment, args.output_dir, args.output_prefix,
                               stand_call_conf=args.emit_quality, max_region_length=args.max_region_len,
                               max_seqs_per_region=args.max_seq_per_region, length_dict=None,
                               parsing_mode="parse", region_list=None,
                               handling_mode=args.handling_mode,
                               job_name=args.slurm_job_name,
                               log_prefix=args.slurm_log_prefix,
                               error_log_prefix=args.slurm_error_log_prefix,
                               max_running_jobs=args.slurm_max_running_jobs,
                               max_running_time=args.slurm_max_running_time,
                               max_memmory_per_cpu=args.slurm_max_memmory_per_cpu,
                               modules_list=args.slurm_modules_list,
                               black_list_scaffold_id_file=args.black_list_scaffold_id_file,
                               gvcf_mode=False)
