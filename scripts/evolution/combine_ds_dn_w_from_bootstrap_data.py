#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import EvolutionRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_dir", action="store", dest="input_dir", required=True,
                    help="Input directory with files containing extracted data")
parser.add_argument("-o", "--output_dir", action="store", dest="output_dir", required=True,
                    help="Output directory")

args = parser.parse_args()

EvolutionRoutines.combine_ds_dn_w_from_bootstrap_data(args.input_dir, args.output_dir, use_node_names_if_possible=True)
