#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Collections.General import SynDict


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input file with top hits")
parser.add_argument("-e", "--header", action="store_true", dest="header",
                    help="Header is present in input file")
parser.add_argument("-o", "--output_file", action="store", dest="output", default="stdout",
                    help="Output file")
parser.add_argument("-k", "--family_column", action="store", dest="fam_col", default=1, type=int,
                    help="Family column position(0-based). Default: 1")
parser.add_argument("-a", "--genes_column", action="store", dest="gen_col", default=0, type=int,
                    help="Genes column position(0-based). Default: 0")

args = parser.parse_args()

hit_dict = SynDict()

hit_dict.read(args.input, header=args.header, allow_repeats_of_key=True,
              key_index=args.fam_col, value_index=args.gen_col)

hit_dict.write(args.output, splited_values=True)
