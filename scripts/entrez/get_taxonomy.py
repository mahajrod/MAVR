#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import NCBIRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input",
                    help="Input file with latin names or ids of taxa (one per line)")
parser.add_argument("-o", "--output", action="store", dest="output",
                    help="Output file")
parser.add_argument("-a", "--email", action="store", dest="email", required=True,
                    help="Email used in Entrez queues")
parser.add_argument("-t", "--input_type", action="store", dest="input_type", default="latin",
                    help="Type of input. Allowed: latin(default), id")

args = parser.parse_args()

NCBIRoutines.get_taxonomy_from_id_file(args.input, args.output, args.email, input_type=args.input_type)

"""
Entrez.email = args.email

taxa_list = read_ids(args.input, header=False)
out_file = open(args.out_prefix, "w")
out_file.write("#species\tlineage\n")
if args.input_type == "latin":
    for taxon in taxa_list:
        print "Handling %s" % taxon
        summary = Entrez.read(Entrez.esearch(db="taxonomy", term=taxon))
        if summary:
            id_list = summary["IdList"]
            for id in id_list:
                print "handling %s" % id
                record = Entrez.read(Entrez.efetch(db="taxonomy", id=id, retmode="xml"))
                out_file.write("%s\t%s\t%s\n" % (taxon, record[0]["Rank"], record[0]["Lineage"]))

elif args.input_type == "id":
    for taxon in taxa_list:
        print "Handling %s" % taxon
        record = Entrez.read(Entrez.efetch(db="taxonomy", id=taxon, retmode="xml"))
        out_file.write("%s\t%s\t%s\n" % (taxon, record[0]["Rank"], record[0]["Lineage"]))


out_file .close()
"""