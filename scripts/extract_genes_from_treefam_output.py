#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--in_file", action="store", dest="in_file",
                    help="File with treefam output")
parser.add_argument("-s", "--species", action="store", dest="species",
                    help="Species to extract")
parser.add_argument("-o", "--out_prefix", action="store", dest="out_prefix",
                    help="Prefix of output  files")
args = parser.parse_args()

out_fd = open("%s.txt" % args.out_prefix, "w")
out_fd.write("#gene_family_id\tgenes(%s)\n" % args.species)
id_fd = open("%s_id.txt" % args.out_prefix, "w")
species_name_length = len(args.species)
with open(args.in_file, "r") as in_fd:
    for line in in_fd:
        tmp_list = line.strip().split("\t")
        family_id = tmp_list[0]
        genes_list = tmp_list[-1].split(",")
        species_genes = []
        for gene in genes_list:
            if gene[-species_name_length:] == args.species:
                species_genes.append(gene[:-species_name_length-1])
                id_fd.write(gene[:-species_name_length-1] + "\n")
        out_fd.write("%s\t%s\n" % (family_id, ','.join(species_genes) if species_genes else "."))

out_fd.close()
id_fd.close()
