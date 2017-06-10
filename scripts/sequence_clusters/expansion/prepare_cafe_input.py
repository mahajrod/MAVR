#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
"""
This scripts assumes that genes are named in following manner: species_geneid
"""

import os
import sys
import argparse

from collections import OrderedDict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()

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
                    help="Input file with hclust output")
parser.add_argument("-c", "--cafe_file", action="store", dest="cafe_file", required=True,
                    help="CAFE file")
parser.add_argument("-s", "--species_set", action="store", dest="species_set", type=lambda s: s.split(","),
                    required=True,
                    help="Comma separated set of species.")

parser.add_argument("-l", "--name_last", action="store_false", dest="name_first", default=True,
                    help="Position of name of species in gene_id")

parser.add_argument("-e", "--name_separator", action="store", dest="name_separator", default="_",
                    help="Separator between species name and gene name. Default: '_'")

args = parser.parse_args()

species_list = sorted(args.species_set)
out_fd = open(args.cafe_file, "w")
out_fd.write("FAMILYDESC\tFAMILY\t%s\n" % ("\t".join(species_list)))


with open(args.input, "r") as in_fd:
    for line in in_fd:
        temp = line.strip().split("\t")
        family_id = temp[0]
        genes = temp[-1][:-1].split(",")

        genes_dict = split_gene_names(genes, name_first=args.name_first, separator=args.name_separator)
        species_genes_dict = OrderedDict()
        out_fd.write("%s\t%s" % (family_id, family_id))
        for species in species_list:
            number_of_genes = len(genes_dict[species]) if species in genes_dict else 0
            # species_genes_dict[species] = number_of_genes
            out_fd.write("\t%i" % number_of_genes)
        out_fd.write("\n")

out_fd.close()