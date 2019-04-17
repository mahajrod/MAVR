#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os
import sys
import argparse
from RouToolPa.Routines.File import check_path, split_filename
from RouToolPa.Collections.General import SynDict



parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_dir", action="store", dest="input", default="./", type=check_path,
                    help="Input directory with fam files")
parser.add_argument("-s", "--species_set", action="store", dest="species_set", type=lambda s: s.split(","),
                    help="Comma separated set of species.")
parser.add_argument("-b", "--use_basename", action="store_true", dest="use_basename",
                    help="Use basenames of files in input directory as species")
parser.add_argument("-u", "--suffix", action="store", dest="suffix", default=".fam",
                    help="Suffix of fam files. Default: .fam")
parser.add_argument("-o", "--output", action="store", dest="output", default="stdout",
                    help="Suffix of fam files")

args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")
species_list = []
suffix_list = []
if args.use_basename:
    for filename in sorted(os.listdir(args.input)):
        dir, basename, ext = split_filename(filename)
        species_list.append(basename)
        suffix_list.append("%s" % ext)
else:
    species_list = sorted(args.species_set)
    suffix_list = [args.suffix for i in range(0, len(species_list))]

out_fd.write("#species\tnumber_of_families\tnumber_of_proteins\n")
for species, suffix in zip(species_list, suffix_list):
    fam_dict = SynDict()
    fam_dict.read("%s%s%s" % (args.input, species, suffix), separator="\t", split_values=True, values_separator=",",
                  key_index=0, value_index=1)
    out_fd.write("%s\t%i\t%i\n" % (species, len(fam_dict), fam_dict.count_all_synonyms()))

if args.output != "stdout":
    out_fd.close()
