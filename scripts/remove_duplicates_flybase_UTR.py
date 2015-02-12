#!/usr/bin/env python
__author__ = 'mahajrod'
import os
import argparse

from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_fasta", action="store", dest="in_fasta",
                    help="Input fasta file")
parser.add_argument("-o", "--out_fasta", action="store", dest="out_fasta",
                    help="Output fasta file")
parser.add_argument("-d", "--index_file", action="store", dest="idx_fasta", default="temp.idx",
                    help="Biopython index of input fasta (will be created is abbsent)")

args = parser.parse_args()

UTR_record_dict = SeqIO.index_db(args.idx_fasta, args.in_fasta, format="fasta")

print("Totally UTRs(with duplicates): %i" % len(UTR_record_dict))

id_list, location_list = [], []
for record_id in UTR_record_dict:
    description = UTR_record_dict[record_id].description
    location_string = description.split(";")[1][1:]
    if id_list:
        if location_list[-1] != location_string:
            id_list.append(record_id)
            location_list.append(location_string)
    else:
        id_list.append(record_id)
        location_list.append(location_string)

print("Totally UTRs(without duplicates): %i" % len(id_list))


def filtered_record_generator(id_list):
    for record_id in id_list:
        yield UTR_record_dict[record_id]

SeqIO.write(filtered_record_generator(id_list), args.out_fasta, "fasta")

os.remove(args.idx_fasta)