#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Evolution import MCMCTree

parser = argparse.ArgumentParser()

parser.add_argument("-s", "--input", action="store", dest="input", required=True,
                    help="Input file alignment")
parser.add_argument("-c", "--ctl_file", action="store", dest="ctl_file", default="mcmctree.ctl",
                    help="File to write configuration. Default: mcmctree.ctl")
parser.add_argument("-t", "--tree_file", action="store", dest="tree",
                    help="File with phylogenetic tree with calibrations")
parser.add_argument("-o", "--output", action="store", dest="output",
                    help="File to write output")
parser.add_argument("--seq_type", action="store", dest="seq_type", default="nucleotides",
                    help="Type of input sequences. Allowed: nucleotides, codons, aminoacids."
                         "Default: nucleotides")
parser.add_argument("--root_age", action="store", dest="root_age",
                    help="Root age, set if it is not present in tree file")
parser.add_argument("--num_of_partitions", action="store", dest="num_of_partitions", type=int,
                    default=1,
                    help="Number of partitions in alignment. Default: 1")
parser.add_argument("--clock_type", action="store", dest="clock_type", default="global",
                    help="Clock type for partitions. Allowed: global, independent, correlated."
                         "Default: global")
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

args = parser.parse_args()

MCMCTree.path = args.dir_path
#MCMCTree.timelog = "mcmctree.timelog"
MCMCTree.run(args.input, args.tree, args.output, args.ctl_file, seed=-1, num_of_partitions=args.num_of_partitions,
             seq_type=args.seq_type, use_data=1, clock=args.clock_type, root_age=args.root_age,
             model=args.model, ncatG=5, alpha_for_gamma_rates_at_sites=args.alpha_for_gamma_rates_at_sites,
             birth=1, death=1, sampling=0.1,
             alpha_gamma_alpha=1, alpha_gamma_beta=1, kappa_gamma_alpha=6, kappa_gamma_beta=2,
             rgene_gamma_alpha=args.rgene_gamma_alpha, rgene_gamma_beta=args.rgene_gamma_beta,
             sigma2_gamma_alpha=1, sigma2_gamma_beta=10,
             remove_ambiguity_sites=args.remove_ambiguity_sites,
             auto_finetune=True, times=0.1, rates=0.1, mixing=0.1, paras=0.1, RateParas=0.1, FossilErr=0.1,
             num_of_burning=args.num_of_burning, sampling_frequency=args.sampling_frequency,
             number_of_samples=args.number_of_samples)




