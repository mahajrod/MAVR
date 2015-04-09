#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os, sys
from collections import OrderedDict

from Routines.File import read_ids

species_list = ["cat",
                "cheetah",
                "lion",
                "tiger",
                "dog",
                "human",
                "mouse",
                "opossum"]

sys.path.append("/Users/mahajrod/genetics/MAVR/")

workdir = "/Users/mahajrod/genetics/Projects/Dobrzhansky/project/gene_families/cafe_run_no_selenocystein/"
data_dir = "/Users/mahajrod/genetics/data/"
treefam_output = workdir + "gene_families_of_eightSpeices"
selenocystein_ids_dict = {}

output_cafe = "cafe_input.cafe"

for species in species_list:
    selenocystein_ids_dict[species] = read_ids(data_dir + species + "/" + species + "_selenocystein_proteins_id.t")

os.chdir(workdir)
selenocystein_families_dict = OrderedDict({})
for species in species_list:
    selenocystein_families_dict[species] = set([])

counts_fd = open("counts.t", "w")
#counts_fd.write = ("#family\tbefore\tafter\tchanges\n")
full_selenocystein_fd = open("full_selonocystein_families", "w")

with open(treefam_output, "r") as in_fd:
    with open(output_cafe, "w") as out_fd:
        header = "FAMILYDESC\tFAMILY\t%s\n" % "\t".join(species_list)
        out_fd.write(header)
        for line in in_fd:
            line_list = line.strip().split()
            family_id = line_list[0]
            number_of_genes = int(line_list[-2])
            genes_list = line_list[-1][:-1].split(",")
            genes_dict = OrderedDict()
            filtered_genes_dict = OrderedDict()
            for species in species_list:
                filtered_genes_dict[species] = []
            for gene_entry in genes_list:
                gene, species = gene_entry.split("_")
                if gene in selenocystein_ids_dict[species]:
                    selenocystein_families_dict[species].add(family_id)
                else:
                    filtered_genes_dict[species].append(gene)

            filtered_gene_counts_list = [len(filtered_genes_dict[species]) for species in species_list]
            number_of_filtered_genes = sum(filtered_gene_counts_list)
            counts_fd.write("%s\t%i\t%i\t%i\t%s\n" % (family_id, number_of_genes, number_of_filtered_genes, number_of_genes - number_of_filtered_genes, "no" if number_of_genes == number_of_filtered_genes else "yes"))
            cafe_string = "%s\t%s\t%s\n" % (family_id, family_id, "\t".join(map(str, filtered_gene_counts_list)))
            if number_of_filtered_genes != 0:
                out_fd.write(cafe_string)
            else:
                full_selenocystein_fd.write(line)

full_selenocystein_fd.close()
counts_fd.close()
if __name__ == "__main__":
    pass