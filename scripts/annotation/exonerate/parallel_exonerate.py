#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Tools.Annotation import Exonerate


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input fasta file with sequences")
parser.add_argument("-a", "--target", action="store", dest="target", required=True,
                    help="File with target sequences")
parser.add_argument("-x", "--annotation", action="store", dest="annotation",
                    help="File with query annotation. Usually cds coordinates in transcript for cdna2genome model")
parser.add_argument("-o", "--output", action="store", dest="output", #required=True,
                    help="Prefix of output files")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads to use")
parser.add_argument("-m", "--model", action="store", dest="model", required=True,
                    help="Model to run")
parser.add_argument("-n", "--number_of_results_to_report", action="store",
                    dest="num_of_results_to_report", default=1, type=int,
                    help="Number of results to report for each input sequence")
parser.add_argument("-u", "--num_of_seq_per_file", action="store", dest="num_of_seq_per_file",
                    type=int, default=None,
                    help="Number of sequences per splited input")
parser.add_argument("-d", "--exonerate_dir", action="store", dest="exonerate_dir", default="",
                    help="Directory with exonerate binary")
parser.add_argument("-e", "--num_of_splited_files", action="store", dest="num_of_splited_files",
                    type=int, default=None,
                    help="Number of splited files")
parser.add_argument("-q", "--softmasked_input", action="store_true", dest="softmasked_input", default=False,
                    help="Input sequences are softmasked")
parser.add_argument("-s", "--softmasked_target", action="store_true", dest="softmasked_target", default=False,
                    help="Target is softmasked")

parser.add_argument("-p", "--handling_mode", action="store", dest="handling_mode", default="local",
                    help="Handling mode. Allowed: local(default), slurm")
parser.add_argument("-j", "--slurm_job_name", action="store", dest="slurm_job_name", default="JOB",
                    help="Slurm job name. Default: JOB")
parser.add_argument("-b", "--slurm_max_jobs", action="store", dest="slurm_max_jobs", default=1000, type=int,
                    help="Slurm max jobs. Default: 1000")
parser.add_argument("-y", "--slurm_log_prefix", action="store", dest="slurm_log_prefix",
                    help="Slurm log prefix. ")
parser.add_argument("-z", "--slurm_cmd_log_file", action="store", dest="slurm_cmd_log_file",
                    help="Slurm cmd logfile")
parser.add_argument("-l", "--slurm_error_log_prefix", action="store", dest="slurm_error_log_prefix",
                    help="Slurm error log prefix")
parser.add_argument("-k", "--max_memory_per_task", action="store", dest="max_memory_per_task", default="5000",
                    help="Maximum memory per task in megabytes. Default: 5000")

parser.add_argument("-r", "--slurm_max_running_time", action="store", dest="slurm_max_running_time", default="100:00:00",
                    help="Slurm max running time in hh:mm:ss format. Default: 100:00:00")

parser.add_argument("-w", "--slurm_modules_list", action="store", dest="slurm_modules_list", default=[],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of modules to load. Set modules for hmmer and python")

parser.add_argument("--splited_fasta_dir", action="store", dest="splited_fasta_dir", default="splited_fasta_dir/",
                    help="Directory to write splited fasta input. Default: splited_fasta_dir/")
parser.add_argument("--splited_output", action="store", dest="splited_output", default="splited_output/",
                    help="Directory to write splited output. Default: splited_output/")
parser.add_argument("--converted_output", action="store", dest="converted_output", default="converted_output/",
                    help="Directory to write converted output. Default: converted_output/")

"""
parser.add_argument("-u", "--num_in_seq_per_file", action="store", dest="num_in_seq_per_file",
                    type=int, default=1000,
                    help="Number of sequences per splited input")
"""
"""
parser.add_argument("-r", "--strand", action="store", dest="strand", default="both",
                    help="Strand to consider. Possible variants: both, forward, backward."
                         "Default: both")

parser.add_argument("-e", "--other_options", action="store", dest="other_options",
                    help="Other augustus options")
parser.add_argument("-c", "--augustus_config_dir", action="store", dest="config_dir",
                    help="Augustus config dir")
parser.add_argument("-p", "--pfam_hmm3", action="store", dest="pfam_db",
                    help="Pfam database in hmm3 format")
parser.add_argument("-w", "--swissprot_blast_db", action="store", dest="swissprot_db",
                    help="Blast database of swissprot")
parser.add_argument("-m", "--masking", action="store", dest="masking",
                    help="Gff of bed file with masking of repeats")
"""
args = parser.parse_args()

if args.num_of_seq_per_file and args.num_of_splited_files:
    raise ValueError("Options -u/--num_of_seq_per_file and -e/--num_of_splited_files can't be set simultaneously")

if (not args.num_of_seq_per_file) and (not args.num_of_splited_files):
    args.num_of_splited_files = 10 * args.threads

Exonerate.threads = args.threads
Exonerate.path = args.exonerate_dir
Exonerate.parallel_alignment(args.input, args.target, args.model, num_of_files=args.num_of_splited_files,
                             num_of_recs_per_file=args.num_of_seq_per_file,
                             show_alignment=True, show_sugar=True,
                             show_cigar=True,
                             show_vulgar=True, show_query_gff=True, show_target_gff=True,
                             store_intermediate_files=True,
                             annotation_file=args.annotation,
                             splited_fasta_dir=args.splited_fasta_dir,
                             splited_result_dir=args.splited_output,
                             number_of_results_to_report=args.num_of_results_to_report,
                             converted_output_dir=args.converted_output,
                             softmasked_target=args.softmasked_target,
                             softmasked_query=args.softmasked_input,
                             cmd_log_file=args.slurm_cmd_log_file,
                             cpus_per_task=1,
                             handling_mode=args.handling_mode,
                             job_name=args.slurm_job_name,
                             log_prefix=args.slurm_log_prefix,
                             error_log_prefix=args.slurm_error_log_prefix,
                             max_jobs=args.slurm_max_jobs,
                             max_running_time=args.slurm_max_running_time,
                             max_memory_per_node=args.max_memory_per_task,
                             modules_list=args.slurm_modules_list,
                             environment_variables_dict=None,
                             length_thresholds=(1000, 2000),
                             memory_thresholds=(6100, 9000, 20000))

