#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

from Bio import Entrez
from Routines.File import read_ids

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input",
                    help="Input file with latin names of taxa (one per line)")
parser.add_argument("-o", "--output_file_prefix", action="store", dest="out_prefix",
                    help="Prefix of output file")
parser.add_argument("-a", "--email", action="store", dest="email", required=True,
                    help="Email used in Entrez queues")

args = parser.parse_args()

Entrez.email = args.email

taxa_list = read_ids(args.input, header=False)
out_file = open(args.out_prefix, "w")
out_file.write("#species\tlineage\n")

for taxon in taxa_list:
    summary = Entrez.read(Entrez.esearch(db="taxonomy", term=taxon))
    if summary:
        id_list = summary["IdList"]
        for id in id_list:
            record = Entrez.read(Entrez.efetch(db="taxonomy", id=id, retmode="xml"))
            out_file.write("%s\t%s\t%s\n" % (taxon, record[0]["Rank"], record[0]["Lineage"]))

out_file .close()