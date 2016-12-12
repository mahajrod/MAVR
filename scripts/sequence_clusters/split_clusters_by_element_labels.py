#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
"""
This scripts assumes that genes are named in following manner: species_geneid
"""

import os
import argparse
from collections import OrderedDict


def split_gene_names(genes_list, label_first=True, separator="_"):
    species_dict = OrderedDict()
    for gene in genes_list:
        gene_name, species = gene.split(separator)

        if label_first:
            gene_name, species = species, gene_name
        if species not in species_dict:
            species_dict[species] = [gene_name]
        else:
            species_dict[species].append(gene_name)

    return species_dict

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input file with clusters")
parser.add_argument("-s", "--labels_set", action="store", dest="species_set", type=lambda s: s.split(","),
                    required=True,
                    help="Comma separated set of element labels.")
parser.add_argument("-d", "--label_clusters_directory", action="store", dest="species_fam_dir", default="species_fam",
                    help="Directory to output separated files with clusters")
parser.add_argument("-l", "--label_last", action="store_true", dest="label_last", default=False,
                    help="Label is located after separator in sequence id. Default: before")
parser.add_argument("-e", "--name_separator", action="store", dest="name_separator", default="@",
                    help="Separator between species name and label. Default: '@'")
parser.add_argument("-a", "--cluster_count_file", action="store", dest="family_count_file",
                    default="cluster_count_file.t",
                    help="File to output cluster counts")
parser.add_argument("-w", "--include_clusters_without_genes", action="store_true", dest="include_without_genes",
                    help="Include clusters without genes in species. Default: false")
parser.add_argument("-c", "--count_sequence", action="store_true", dest="count_genes",
                    help="Count sequences in clusters")
parser.add_argument("-r", "--header", action="store_true", dest="header",
                    help="Header is present in input file with clusters")
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

        genes = split_gene_names(genes, label_first=not args.label_last, separator=args.name_separator)

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


