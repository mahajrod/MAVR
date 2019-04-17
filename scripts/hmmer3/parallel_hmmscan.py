#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.HMMER import HMMER3
from RouToolPa.Routines.File import check_path


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
parser.add_argument("-d", "--output_dir", action="store", dest="output_dir",
                    default="./", type=check_path,
                    help="Directory to write output. Default: current directory")
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

parser.add_argument("-u", "--slurm_max_memmory_per_cpu", action="store", dest="slurm_max_memmory_per_cpu",
                    default=4000, type=int,
                    help="Slurm maximum memmory per cpu in megabytes. Default: 4000")
parser.add_argument("-q", "--slurm_modules_list", action="store", dest="slurm_modules_list", default=[],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of modules to load. Set modules for hmmer and python")

parser.add_argument("-r", "--MAVR_script_dir", action="store", dest="MAVR_script_dir", default="",
                    help="MAVR scrip dir. Default: ''")

parser.add_argument("-g", "--parsing_mode", action="store", dest="parsing_mode", default="index_db",
                    help="Parsing mode for hmmer hits file. Allowed: parse, index, index_db(default)")

args = parser.parse_args()


HMMER3.threads = 1
HMMER3.path = args.path
HMMER3.parallel_hmmscan(args.input, args.input_seq, args.output_prefix, args.output_dir,
                        num_of_seqs_per_scan=None,
                        threads=args.threads,
                        combine_output_to_single_file=args.combine_output,
                        dont_output_alignments=args.no_alignment,
                        handling_mode=args.handling_mode,
                        job_name=args.slurm_job_name,
                        log_prefix=args.slurm_log_prefix,
                        error_log_prefix=args.slurm_error_log_prefix,
                        max_running_jobs=args.slurm_max_running_jobs,
                        max_running_time=args.slurm_max_running_time,
                        max_memmory_per_cpu=args.slurm_max_memmory_per_cpu,
                        MAVR_scripts_dir=args.MAVR_script_dir,
                        hmm_hit_parsing_mode=args.parsing_mode,
                        modules_list=args.slurm_modules_list,
                        extract_top_hits=False,
                        get_clusters_from_top_hits=False
                        )

