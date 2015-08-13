#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os
import sys
import argparse

from collections import OrderedDict

from Routines.File import check_path
from CustomCollections.GeneralCollections import SynDict, TwoLvlDict, IdList

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_dir", action="store", dest="input", required=True, type=check_path,
                    help="Input directory with fam files")
parser.add_argument("-s", "--species_set", action="store", dest="species_set", type=lambda s: s.split(","),
                    required=True,
                    help="Comma separated set of species.")
parser.add_argument("-u", "--suffix", action="store", dest="suffix", required=True, default=".fam",
                    help="Suffix of fam files. Default: .fam")
parser.add_argument("-o", "--output", action="store", dest="output", default="stdout",
                    help="Suffix of fam files")

args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")
species_list = sorted(args.species_set)

out_fd.write("#species\tnumber_of_families\tnumber_of_proteins\n")
for species in species_list:
    fam_dict = SynDict()
    fam_dict.read("%s%s" % (species, args.suffix), separator="\t", split_values=False, values_separator=",",
                  key_index=0, value_index=1)
    out_fd.write("%s\t%i\t%i\n" % (species, len(fam_dict), fam_dict.count_all_synonyms()))

if args.output != "stdout":
    out_fd.close()
