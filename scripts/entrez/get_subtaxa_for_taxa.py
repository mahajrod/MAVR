#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import NCBIRoutines
from RouToolPa.Collections.General import IdList


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input_file",
                    help="Input file with latin names or taxa (one per line)")
parser.add_argument("-t", "--taxa_list", action="store", dest="taxa_list",
                    type=lambda s: s.split(","),
                    help="Comma-separated list of taxa latin names")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output file")
parser.add_argument("-a", "--email", action="store", dest="email", required=True,
                    help="Email used in Entrez queues")
parser.add_argument("-s", "--subtaxa_rank_list", action="store", dest="subtaxa_rank_list",
                    default=['species'], type=lambda s: s.split(","),
                    help="Comma-separated list of subtaxa ranksl to extract. Default: ['species']")
parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                    help="Verbose log")

args = parser.parse_args()

if args.input_file and args.taxa_list:
    raise ValueError("ERROR! Both input file and taxa list were set")
elif args.input_file:
    taxa_list = IdList(filename=args.input_file)
elif args.taxa_list:
    taxa_list = args.taxa_list
else:
    raise ValueError("ERROR! Neither input file nor taxa list was set")

NCBIRoutines.get_subtaxa_for_taxa(taxa_list,
                                  args.email,
                                  args.output_prefix,
                                  subtaxa_rank=args.subtaxa_rank_list,
                                  verbose=args.verbose)
