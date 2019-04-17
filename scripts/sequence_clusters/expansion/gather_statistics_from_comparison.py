#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from collections import OrderedDict
from RouToolPa.Collections.General import TwoLvlDict
from RouToolPa.Routines.File import check_path



parser = argparse.ArgumentParser()

parser.add_argument("-s", "--species_list", action="store", dest="species_list", type=lambda s: s.split(","),
                    required=True,
                    help="Comma-separated list of species")
parser.add_argument("-d", "--species_dir", action="store", dest="species_dir", default="./", type=check_path,
                    help="Directory with per species statistics")
parser.add_argument("-o", "--output_file", action="store", dest="output", default="stdout",
                    help="Output file. Default: stdout")

args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

species_stat_dict = TwoLvlDict()

for species in args.species_list:
    with open("%s%s/stat.t" % (args.species_dir, species), "r") as stat_fd:
        statistics = map(lambda s: s.strip().split("\t"), stat_fd.readlines())
    species_stat_dict[species] = OrderedDict(statistics)

species_stat_dict.write(out_fd)
if args.output != "stdout":
    out_fd.close()
