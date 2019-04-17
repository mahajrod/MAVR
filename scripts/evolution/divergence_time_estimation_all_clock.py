#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Evolution import MCMCTree

parser = argparse.ArgumentParser()

parser.add_argument("-s", "--input", action="store", dest="input", required=True,
                    help="Input file alignment")
parser.add_argument("-t", "--tree_file", action="store", dest="tree",
                    help="File with phylogenetic tree with calibrations")
parser.add_argument("-o", "--outdir", action="store", dest="outdir", default="./",
                    help="Output directory. Default: ./")
parser.add_argument("-p", "--output_prefix", action="store", dest="output_prefix", default="divtime",
                    help="Prefix of output files. Default: divtime")
parser.add_argument("--seq_type", action="store", dest="seq_type", default="nucleotides",
                    help="Type of input sequences. Allowed: nucleotides, codons, aminoacids."
                         "Default: nucleotides")
parser.add_argument("--root_age", action="store", dest="root_age",
                    help="Root age, set if it is not present in tree file")
parser.add_argument("--num_of_partitions", action="store", dest="num_of_partitions", type=int,
                    default=1,
                    help="Number of partitions in alignment. Default: 1")
parser.add_argument("-m", "--model", action="store", dest="model", required=True,
                    help="Model of substitutions to use. "
                         "Allowed: JC69, K80, F81, F84, HKY85") #, TC92, TN93, GTR, UNREST, REVu, UNSETu")
parser.add_argument("--alpha_for_gamma_rates_at_sites", action="store", dest="alpha_for_gamma_rates_at_sites",
                    type=float, default=0.5,
                    help="Alpha for gamma rates at sites. If alpha != 0, the program will assume a gamma-rates "
                         "model, while alpha = 0 means that the model of one rate for all sites will be used. ")
parser.add_argument("--num_of_burning", action="store", dest="num_of_burning", type=int,
                    default=2000,
                    help="Number of generation to be counted as burning phase. Default: 2000")
parser.add_argument("--sampling_frequency", action="store", dest="sampling_frequency", type=int,
                    default=2,
                    help="Sampling frequency(sample every NUM generations). Default: 2")
parser.add_argument("--number_of_samples", action="store", dest="number_of_samples", type=int,
                    default=20000,
                    help="Number of samples retrieved from generations. Default: 20000")
parser.add_argument("--remove_ambiguity_sites", action="store_true", dest="remove_ambiguity_sites",
                    default=False,
                    help="Remove ambiguity sites")
parser.add_argument("--rgene_gamma_alpha", action="store", dest="rgene_gamma_alpha", type=float, default=2.0,
                    help="Specifies the shape(alpha) parameter in the gamma prior for the overall rate parameter mu"
                         "Default: 2")
parser.add_argument("--rgene_gamma_beta", action="store", dest="rgene_gamma_beta", type=float, default=2.0,
                    help="Specifies the scale(beta) parameter in the gamma prior for the overall rate parameter mu"
                         "Default: 2")
parser.add_argument("-d", "--directory_with_mcmctree_bin", action="store", dest="dir_path",
                    help="Path to directory with mcmctree binary")

parser.add_argument("-q", "--handling_mode", action="store", dest="handling_mode", default="local",
                    help="Handling mode. Allowed: local(default), slurm")
parser.add_argument("-j", "--slurm_job_name", action="store", dest="slurm_job_name", default="JOB",
                    help="Slurm job name. Default: JOB")
#parser.add_argument("-m", "--slurm_max_jobs", action="store", dest="slurm_max_jobs", default=1000, type=int,
#                    help="Slurm max jobs. Default: 1000")
"""
parser.add_argument("-y", "--slurm_log_prefix", action="store", dest="slurm_log_prefix",
                    help="Slurm log prefix. ")
parser.add_argument("-z", "--slurm_cmd_log_file", action="store", dest="slurm_cmd_log_file",
                    help="Slurm cmd logfile")
parser.add_argument("-l", "--slurm_error_log_prefix", action="store", dest="slurm_error_log_prefix",
                    help="Slurm error log prefix")
"""
parser.add_argument("-e", "--max_memory_per_task", action="store", dest="max_memory_per_task", default="3000",
                    help="Maximum memory per task in megabytes. Default: 3000")

parser.add_argument("-a", "--slurm_max_running_time", action="store", dest="slurm_max_running_time", default="100:00:00",
                    help="Slurm max running time in hh:mm:ss format. Default: 100:00:00")

parser.add_argument("-w", "--slurm_modules_list", action="store", dest="slurm_modules_list", default=[],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of modules to load. Set modules for hmmer and python")

args = parser.parse_args()

MCMCTree.path = args.dir_path
MCMCTree.run_all_clocks(args.input, args.tree, args.outdir, output_prefix=args.output_prefix, seed=-1,
                        num_of_partitions=args.num_of_partitions,
                        seq_type=args.seq_type, use_data=1, root_age=args.root_age,
                        model=args.model, ncatG=5,
                        alpha_for_gamma_rates_at_sites=args.alpha_for_gamma_rates_at_sites,
                        birth=1, death=1, sampling=0.1,
                        alpha_gamma_alpha=1.0, alpha_gamma_beta=1.0, kappa_gamma_alpha=6.0, kappa_gamma_beta=2.0,
                        rgene_gamma_alpha=args.rgene_gamma_alpha, rgene_gamma_beta=args.rgene_gamma_beta,
                        sigma2_gamma_alpha=1.0, sigma2_gamma_beta=10.0,
                        remove_ambiguity_sites=args.remove_ambiguity_sites,
                        auto_finetune=True, times=0.1, rates=0.1, mixing=0.1, paras=0.1, RateParas=0.1, FossilErr=0.1,
                        num_of_burning=args.num_of_burning,
                        sampling_frequency=args.sampling_frequency,
                        number_of_samples=args.number_of_samples,
                        cpus_per_task=1,
                        handling_mode=args.handling_mode,
                        job_name=args.slurm_job_name,
                        max_jobs=None,
                        max_running_time=args.slurm_max_running_time,
                        max_memory_per_node=args.max_memory_per_task,
                        max_memmory_per_cpu=None,
                        modules_list=args.slurm_modules_list,
                        environment_variables_dict=None,
                        duplicate_to_stdout=True)

"""
MCMCTree.run(args.input, args.tree, args.output, args.ctl_file, seed=-1, num_of_partitions=args.num_of_partitions,
             seq_type=args.seq_type, use_data=1, clock=args.clock_type, root_age=args.root_age,
             model=args.model, ncatG=5, alpha_for_gamma_rates_at_sites=args.alpha_for_gamma_rates_at_sites,
             birth=1, death=1, sampling=0.1,
             alpha_gamma_alpha=1, alpha_gamma_beta=1, kappa_gamma_alpha=6, kappa_gamma_beta=2,
             rgene_gamma_alpha=args.rgene_gamma_alpha, rgene_gamma_beta=args.rgene_gamma_beta,
             sigma2_gamma_alpha=1, sigma2_gamma_beta=10,
             remove_ambiguity_sites=args.remove_ambiguity_sites,
             auto_finetune=True, times=0.1, rates=0.1, mixing=0.1, paras=0.1, RateParas=0.1, FossilErr=0.1,
             num_of_burning=args.num_of_burning, sampling_frequency=args.number_of_samples,
             number_of_samples=args.number_of_samples)
"""""



