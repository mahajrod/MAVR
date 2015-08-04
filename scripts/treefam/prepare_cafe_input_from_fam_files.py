#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
"""
This scripts assumes that genes are named in following manner: species_geneid
"""

import os
import sys
import argparse

from collections import OrderedDict
from CustomCollections.GeneralCollections import SynDict, TwoLvlDict

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_dir", action="store", dest="input", required=True,
                    help="Input directory with fam files. Names of fam files should"
                         "start with species name and end with suffix set by -u option")
parser.add_argument("-c", "--cafe_file", action="store", dest="cafe_file", required=True,
                    help="CAFE file")
parser.add_argument("-s", "--species_set", action="store", dest="species_set", type=lambda s: s.split(","),
                    required=True,
                    help="Comma separated set of species.")
parser.add_argument("-u", "--suffix", action="store", dest="suffix", required=True,
                    help="Suffix of fam files")

args = parser.parse_args()

species_list = sorted(args.species_set)
out_fd = open(args.cafe_file, "w")
out_fd.write("FAMILYDESC\tFAMILY\t%s\n" % ("\t".join(species_list)))

fam_count_dict = TwoLvlDict()

for species in args.species_set:
    species_fam = SynDict()
    species_fam.read("%s%s" % (species, args.suffix), split_values=True,
                     values_separator=",", separator="\t")
    fam_count_dict[species] = species_fam.count_synonyms()

for family in fam_count_dict.sl_keys():
    number_str = "\t".join(map(str,
                               [fam_count_dict[species][family] if family in fam_count_dict[species] else 0
                                for species in species_list]))
    out_fd.write("%s\t%s\t%s\n" % (family, family, number_str))

