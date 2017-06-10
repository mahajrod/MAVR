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

parser.add_argument("-i", "--input_file", action="store", dest="input", default="stdin",
                    help="Input file. Default: stdin")
parser.add_argument("-o", "--output_file", action="store", dest="output", default="stat.svg",
                    help="Output file with pictures. Default: stat.svg")
parser.add_argument("-n", "--number_of_species", action="store", dest="number_of_species", type=int,
                    help="Number of species in analysis.")
parser.add_argument("-s", "--species_set", action="store", dest="species_set",
                    help="Comma-separated list of species.")
parser.add_argument("-l", "--name_last", action="store_false", dest="name_first", default=True,
                    help="Position of name of species in gene_id")
parser.add_argument("-e", "--name_separator", action="store", dest="name_separator", default="_",
                    help="Separator between species name and gene name. Default: '_'")
args = parser.parse_args()

args.species_set = set(args.species_set.split(","))

in_fd = sys.stdin if args.input == "stdin" else open(args.input, "r")


#n_x_bins = np.linspace(1, 40, 40 - int(size) + 1)
#n_y_bins = np.linspace(1, 40, 40 - int(size) + 1)
number_of_genes_in_family_list = []
families_list = []
species_in_family_list = []
species_families_dict = OrderedDict()
for species in args.species_set:
    species_families_dict[species] = []
for line in in_fd:
    gene_list = line.strip().split("\t")[6][:-1].split(",")
    number_of_genes_in_family_list.append(int(line.strip().split("\t")[5]))
    species_dict = split_gene_names(gene_list, name_first=args.name_first, separator=args.name_separator)
    families_list.append(species_dict)
    species_in_family_list.append(len(species_dict))
    for species in args.species_set:
        if species in species_dict:
            species_families_dict[species].append(len(species_dict[species]))
max_number_family = max(number_of_genes_in_family_list)

figure = plt.figure(1, figsize=(10, 10))
ax1 = figure.add_subplot(2, 2, 1)
ax1_bins = np.arange(1, 101, 1)
ax1_bins = np.append(ax1_bins, max(number_of_genes_in_family_list))
hist, bin_edges = np.histogram(number_of_genes_in_family_list, bins=ax1_bins)
plt.bar(bin_edges[:-1], hist, width=1)
plt.xticks((1, 20, 40, 60, 80, 100), ('1', '20', '40', '60', '80', '100+'))
plt.xlim(xmin=1)

plt.xlabel("N of genes in cluster")
plt.ylabel("N of clusters")
plt.xlim(xmax=101)
ax2 = plt.subplot(2, 2, 2)
#print(species_in_family_list)
bins = np.arange(1, args.number_of_species + 2, 1)
hist_out = plt.hist(species_in_family_list, bins=bins, align='left', color='green')

plt.xlabel("Species with genes in cluster")
plt.ylabel("N of clusters")
plt.xticks([i for i in range(1, args.number_of_species + 1)], [str(i) for i in range(1, args.number_of_species + 1)])
plt.xlim(xmin=0.5, xmax=args.number_of_species + 0.5)


ax3 = plt.subplot(2, 1, 2)

hist_dict = OrderedDict()
plot_list = []
ax3_bins = np.arange(1, 26, 1)
ax3_bins = np.append(ax3_bins, max(number_of_genes_in_family_list))

color_list = ["blue",
              "green",
              "red",
              "cyan",
              "magenta",
              "orange",
              "black",
              "teal",
              "darkred",
              "yellowgreen"]
index = 0
for species in species_families_dict:
    hist_dict[species] = np.histogram(species_families_dict[species], bins=ax3_bins)
    plt.plot(hist_dict[species][1][:-1], hist_dict[species][0], label=species, color=color_list[index])
    index += 1

plt.legend()
plt.xlim(xmin=1)
plt.xticks((1, 5, 10, 15, 20, 25), ('1', '5', '10', '15', '20', '25+'))
plt.xlabel("N of genes in cluster")
plt.ylabel("N of clusters")
plt.subplots_adjust(hspace=0.25, wspace=0.25)
plt.savefig(args.output)

if args.input != "stdin":
    in_fd.close()