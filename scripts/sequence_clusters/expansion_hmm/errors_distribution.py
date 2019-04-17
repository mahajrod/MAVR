#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.Collections.General import TwoLvlDict




def is_assembled(families):
    fam_set = list(set(families))
    if len(fam_set) > 1 or "_" in fam_set[0]:
        return False
    return True

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="File with families assembly information")
parser.add_argument("-e", "--header", action="store_true", dest="header",
                    help="Header is present in input file")
parser.add_argument("-o", "--output_file", action="store", dest="output", default="stdout",
                    help="Output file. Default: stdout")

args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

species_syn_dict = TwoLvlDict()
out_fd.write("#family\tspecies_with_family\tspecies_with_errors\tspecies_with_correct_fam\terror_ratio\n")
with open(args.input, "r") as in_fd:
    if args.header:
        in_fd.readline()
    for line in in_fd:
        species_with_errors = 0
        species_with_fam = 0
        tmp = line.strip().split("\t")
        family_name = tmp[0]
        for fam in tmp[1:]:
            if fam != ".":
                species_with_fam += 1
            if "_" in fam:
                species_with_errors += 1
        species_with_correct_fam = species_with_fam - species_with_errors
        error_ratio = float(species_with_errors)/float(species_with_fam)
        out_fd.write("%s\t%i\t%i\t%i\t%.4f\n" %
                     (family_name, species_with_fam, species_with_errors, species_with_correct_fam, error_ratio))

if args.output != "stdout":
    out_fd.close()
