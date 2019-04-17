#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import PhylogeneticsRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input alignment")
parser.add_argument("-o", "--output_dir", action="store", dest="output_dir", required=True,
                    help="Output directory")
parser.add_argument("-p", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Output prefix")
parser.add_argument("-n", "--replicate_number", action="store", dest="replicate_number", type=int,
                    help="Number of bootstrap replicates")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of alignment. Default: fasta")


args = parser.parse_args()

PhylogeneticsRoutines.bootstrap_alignment(args.input, args.output_dir, args.output_prefix,
                                          args.replicate_number, format=args.format)
