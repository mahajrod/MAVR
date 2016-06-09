#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse

from Tools.Evolution import Codeml

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_dir", action="store", dest="input_dir", required=True,
                    help="Input directory with codon alignments")
parser.add_argument("-r", "--tree_file", action="store", dest="tree", required=True,
                    help="File with phylogenetic tree in Newick format")
parser.add_argument("-o", "--output_dir", action="store", dest="output_dir", required=True,
                    help="Output directory")
parser.add_argument("--seq_type", action="store", dest="seq_type", default="codons",
                    help="Type of input sequences. Allowed: codons, aminoacids, translate."
                         "Default: codons")
parser.add_argument("-c", "--codon_frequency", action="store", dest="codon_frequency", default="F3X4",
                    help="Codon frequency. Allowed: equal, F1X4, F3X4, table. Default: F3X4")
parser.add_argument("-g", "--genetic_code", action="store", dest="genetic_code", default=0,
                    help="Genetic code. Allowed: 0-10. Default: 0 - mammalian.")
parser.add_argument("-m", "--model", action="store", dest="model", type=int, default=1,
                    help="Model to use for dN/dS calculations. Allowed: 0(same for all branches), "
                         "1(independent for each branch),"
                         "2(2 or more for branches). Default - 1")
parser.add_argument("-d", "--dont_clean_data", action="store_false", dest="clean_data", default=True,
                    help="Don't clean input data")
parser.add_argument("-s", "--small_difference", action="store", dest="small_difference",
                    type=float, default=0.00001,
                    help="Maximum difference to stop simulations. Default - 0.00001")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int,
                    default=1,
                    help="Number of threads to use")
parser.add_argument("-p", "--path", action="store", dest="path",
                    help="Path to directory with PAML binaries")
args = parser.parse_args()

Codeml.threads = args.threads
print Codeml.threads
Codeml.path = args.path
Codeml.parallel_codeml(args.input_dir, args.tree, args.output_dir, seq_type=args.seq_type,
                       codon_frequency=args.codon_frequency, noisy=3, verbose="concise", runmode=0, clock=0,
                       aminoacid_distance=None, model=args.model, nssites=0,
                       genetic_code=args.genetic_code, fix_kappa=False, kappa=5, fix_omega=False, omega=0.2, getSE=0,
                       RateAncestor=0, small_difference=args.small_difference, clean_data=args.clean_data, method=0)




