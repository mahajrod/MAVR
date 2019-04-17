#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import EvolutionRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with bootstrap samples")
parser.add_argument("-o", "--output_dir", action="store", dest="output_dir", required=True,
                    help="Output directory")
parser.add_argument("-p", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Output prefix")
parser.add_argument("-e", "--out_ext", action="store", dest="out_ext", default="phy",
                    help="Extension of output files. Default: phy")

args = parser.parse_args()

EvolutionRoutines.split_paml_bootstrap_samples(args.input, args.output_dir, args.output_prefix,
                                               output_extension=args.out_ext)
