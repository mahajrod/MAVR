#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import argparse

from collections import OrderedDict
from CustomCollections.GeneralCollections import TwoLvlDict

from Routines.File import check_path
from Routines.File import read_synonyms_dict

parser = argparse.ArgumentParser()

parser.add_argument("-s", "--species_list", action="store", dest="species_list", type=lambda s: s.split(","),
                    required=True,
                    help="Comma-separated list of species")
parser.add_argument("-d", "--species_dir", action="store", dest="species_dir", default="./", type=check_path,
                    help="Comma-separated list of species")
parser.add_argument("-o", "--output_file", action="store", dest="output", default="stdout",
                    help="Output file. Default: stdout")

args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

species_syn_dict = TwoLvlDict()

for species in args.species_list:
    species_syn_dict[species] = read_synonyms_dict("%s%s/synonym.t" % (args.species_dir, species))

species_syn_dict.write("families_all_species.t", absent_symbol=".")
if args.output != "output":
    out_fd.close()
