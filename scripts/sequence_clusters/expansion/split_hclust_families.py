#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
"""
This scripts assumes that genes are named in following manner: species_geneid
"""

import os
import argparse

from collections import OrderedDict


def split_gene_names(genes_list, name_first=True, separator="_"):
    species_dict = OrderedDict()
    for gene in genes_list:
        gene_name, species = gene.split(separator)

        if name_first:
            gene_name, species = species, gene_name
        if species not in species_dict:
            species_dict[species] = [gene_name]
        else:
            species_dict[species].append(gene_name)

    return species_dict

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input file with hclust families")
parser.add_argument("-s", "--species_set", action="store", dest="species_set", type=lambda s: s.split(","),
                    required=True,
                    help="Comma separated set of species.")
parser.add_argument("-d", "--species_families_directory", action="store", dest="species_fam_dir", default="species_fam",
                    help="Directory to output files with species families")
parser.add_argument("-l", "--name_last", action="store_false", dest="name_first", default=True,
                    help="Position of name of species in gene_id")
parser.add_argument("-e", "--name_separator", action="store", dest="name_separator", default="_",
                    help="Separator between species name and gene name. Default: '_'")
parser.add_argument("-a", "--family_count_file", action="store", dest="family_count_file",
                    default="family_count_file.t",
                    help="File to output family counts")
parser.add_argument("-w", "--include_families_without_genes", action="store_true", dest="include_without_genes",
                    help="Include famalies without genes in species")
parser.add_argument("-c", "--count_genes", action="store_true", dest="count_genes",
                    help="Count genes in families")
parser.add_argument("-r", "--header", action="store_true", dest="header",
                    help="Header is present in input file with families")
args = parser.parse_args()

if len(args.name_separator) > 1:
    raise ValueError("Length of name separator is not equal 1.")

try:
    os.mkdir(args.species_fam_dir)
except OSError:
    pass

species_fd_dict = OrderedDict()
family_count_dict = OrderedDict()
for species in args.species_set:
    species_fd_dict[species] = open("%s/%s.fam" % (args.species_fam_dir, species), "w")
    family_count_dict[species] = 0

all_fam_fd = open("all_species.fam", "w")
with open(args.input, "r") as in_fd:
    if args.header:
        in_fd.readline()
    for line in in_fd:
        temp = line.strip().split("\t")
        family_id = temp[0]
        all_fam_fd.write("%s\t%s\n" % (family_id, temp[-1]))
        genes = temp[-1].split(",")

        genes = split_gene_names(genes, name_first=args.name_first, separator=args.name_separator)

        if args.include_without_genes:
            for species in args.species_set:
                if species not in genes:
                    genes[species] = ["."]

        for species in genes:
            if args.count_genes:
                species_fd_dict[species].write("%s\t%i\t%s\n" % (family_id, len(genes[species]),
                                                                 ",".join(genes[species])))
            else:
                species_fd_dict[species].write("%s\t%s\n" % (family_id, ",".join(genes[species])))
            family_count_dict[species] += 1
for species in species_fd_dict:
    species_fd_dict[species].close()
all_fam_fd.close()

with open(args.family_count_file, "w") as count_fd:
    count_fd.write("#species\tnumber_of_families\n")
    for species in family_count_dict:
        count_fd.write("%s\t%i\n" % (species, family_count_dict[species]))


