#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse

from Tools.HMMER import HMMER3
from Routines.File import check_path

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_hmm", action="store", dest="input",
                    help="Input hmm3 file")
parser.add_argument("-s", "--input_seq", action="store", dest="input_seq",
                    help="Input file with sequences")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    help="Prefix of output files")
parser.add_argument("-c", "--combine_output", action="store_true", dest="combine_output",
                    help="Combine output files to single")
parser.add_argument("--no_ali", action="store_true", dest="no_alignment",
                    help="Dont save alignments to minimize output")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads")
parser.add_argument("-d", "--hmmscan_output_dir", action="store", dest="hmmscan_output_dir",
                    default="hmmscan_output_dir/", type=check_path,
                    help="Directory to write intermediate(splited) output")

parser.add_argument("--tblout_dir", action="store", dest="tblout_dir",
                    default="tblout_dir", type=check_path,
                    help="Directory to write intermediate(splited) parseable table of per-sequence hits")
parser.add_argument("--domtblout_dir", action="store", dest="domtblout_dir",
                    default="domtblout_dir", type=check_path,
                    help="Directory to write intermediate(splited) parseable table of per-domain hits")
parser.add_argument("--pfamtblout_dir", action="store", dest="pfamtblout_dir",
                    default="pfamtblout_dir", type=check_path,
                    help="Directory to write intermediate(splited) table of hits and domains to file, in Pfam format ")

#parser.add_argument("--tblout", action="store", dest="tblout",
#                    help="File to save parseable table of per-sequence hits")
#parser.add_argument("--domtblout", action="store", dest="domtblout",
#                    help="File to save parseable table of per-domain hits")
#parser.add_argument("--pfamtblout", action="store", dest="pfamtblout",
#                    help="File to save table of hits and domains to file, in Pfam format ")
parser.add_argument("--hmmer_dir", action="store", dest="path", default="",
                    help="Path to directory with hmmer3.1 binaries")
parser.add_argument("-m", "--handling_mode", action="store", dest="handling_mode", default="local",
                    help="Handling mode. Allowed: local(default), slurm")
parser.add_argument("-j", "--slurm_job_name", action="store", dest="slurm_job_name", default="JOB",
                    help="Slurm job name. Default: JOB")

parser.add_argument("-l", "--slurm_log_prefix", action="store", dest="slurm_log_prefix",
                    help="Slurm log prefix. ")
parser.add_argument("-e", "--slurm_error_log_prefix", action="store", dest="slurm_error_log_prefix",
                    help="Slurm error log prefix")
parser.add_argument("-x", "--slurm_max_running_jobs", action="store", dest="slurm_max_running_jobs",
                    default=300, type=int,
                    help="Slurm max running jobs. Default: 300")
parser.add_argument("-a", "--slurm_max_running_time", action="store", dest="slurm_max_running_time", default="100:00:00",
                    help="Slurm max running time in hh:mm:ss format. Default: 100:00:00")
parser.add_argument("-j", "--slurm_job_name", action="store", dest="slurm_job_name", default="JOB",
                    help="Slurm job name. Default: JOB")
parser.add_argument("-u", "--slurm_max_memmory_per_cpu", action="store", dest="slurm_max_memmory_per_cpu",
                    default=4000, type=int,
                    help="Slurm maximum memmory per cpu in megabytes. Default: 4000")

args = parser.parse_args()


HMMER3.threads = 1
HMMER3.path = args.path
HMMER3.parallel_hmmscan(args.input, args.input_seq, args.output_prefix, "./", num_of_seqs_per_scan=None, split_dir="splited_fasta",
                        splited_output_dir=args.hmmscan_output_dir, threads=args.threads,
                        combine_output_to_single_file=args.combine_output, dont_output_alignments=args.no_alignment,
                        splited_tblout_dir=args.tblout_dir, splited_domtblout_dir=args.domtblout_dir,
                        splited_pfamtblout_dir=args.pfamtblout_dir,
                        handling_mode="local",
                        job_name=args.slurm_job_name,
                        log_prefix=args.slurm_log_prefix,
                        error_log_prefix=args.slurm_error_log_prefix,
                        job_array_script_file=args.slurm_job_array_script_file,
                         #task_index_list=None,
                         #start_task_index=None,
                         #end_task_index=None,
                        max_running_jobs=args.slurm_max_running_jobs,
                        max_running_time=args.slurm_max_running_time,
                        max_memmory_per_cpu=args.slurm_max_memmory_per_cpu
                        )

