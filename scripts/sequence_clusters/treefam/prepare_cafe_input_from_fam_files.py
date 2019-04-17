#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
"""
This scripts assumes that genes are named in following manner: species_geneid
"""

import argparse
from collections import OrderedDict
from RouToolPa.Routines import FileRoutines
from RouToolPa.Collections.General import SynDict, TwoLvlDict, IdList


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
parser.add_argument("-b", "--black_list_file", action="store", dest="black_list_file",
                    help="File with family ids from black list")
parser.add_argument("-w", "--white_list_file", action="store", dest="white_list_file",
                    help="File with family ids from white list")
parser.add_argument("-m", "--min_species_number", action="store", dest="min_species_number", default=1, type=int,
                    help="Minimum number of species with family to retain family. Default: 1")
parser.add_argument("-f", "--filtered_families_directory", action="store", dest="filtered_family_dir",
                    default="filtered_fam", type=FileRoutines.check_path,
                    help="Directory to write filtered_families")
args = parser.parse_args()

FileRoutines.safe_mkdir(args.filtered_family_dir)
species_list = sorted(args.species_set)
if args.white_list_file and args.black_list_file:
    raise ValueError("Black list and white list cant be set simultaneously")

black_list = IdList()
white_list = IdList()
if args.black_list_file:
    black_list.read(args.black_list_file)
if args.white_list_file:
    white_list.read(args.white_list_file)
out_fd = open(args.cafe_file, "w")
filtered_fd = open("%sfiltered_families.cafe" % args.filtered_family_dir, "w")
out_fd.write("FAMILYDESC\tFAMILY\t%s\n" % ("\t".join(species_list)))
filtered_fd.write("FAMILYDESC\tFAMILY\t%s\n" % ("\t".join(species_list)))
species_filtered_fd_list = OrderedDict()
fam_count_dict = TwoLvlDict()
species_family_dict = TwoLvlDict()
for species in args.species_set:
    species_family_dict[species] = SynDict()
    species_family_dict[species].read("%s%s%s" % (FileRoutines.check_path(args.input), species, args.suffix), split_values=True,
                                      values_separator=",", separator="\t")
    #print species_family_dict[species]
    fam_count_dict[species] = species_family_dict[species].count_synonyms()
    #print fam_count_dict[species]
    species_filtered_fd_list[species] = open("%s%s.fam" % (args.filtered_family_dir, species), "w")

for family in fam_count_dict.sl_keys():
    genes_number_list = []
    number_of_species = 0
    for species in species_list:
        genes_number_list.append(fam_count_dict[species][family] if family in fam_count_dict[species] else 0)
        number_of_species += 1 if family in fam_count_dict[species] else 0

    number_str = "\t".join(map(str,
                               genes_number_list))

    if (black_list and (family in black_list)) or (white_list and (family not in white_list)) or (number_of_species < args.min_species_number):
        filtered_fd.write("%s\t%s\t%s\n" % (family, family, number_str))
        for species in species_list:
            if family in species_family_dict[species]:
                species_filtered_fd_list[species].write("%s\t%s\n" % (family, ",".join(species_family_dict[species][family])))
        continue

    out_fd.write("%s\t%s\t%s\n" % (family, family, number_str))

filtered_fd.close()
out_fd.close()
for species in species_filtered_fd_list:
    species_filtered_fd_list[species].close()
