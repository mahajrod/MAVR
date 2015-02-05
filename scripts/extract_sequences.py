#!/usr/bin/env python2
__author__ = 'mahajrod'
import argparse
import os

from Bio import SeqIO
from BCBio import GFF

parser = argparse.ArgumentParser()

parser.add_argument("-g", "--in_gff", action="store", dest="in_gff",
                    help="input gff file")
parser.add_argument("-i", "--in_fasta", action="store", dest="in_fasta",
                    help="input fasta file")
parser.add_argument("-o", "--out_fasta", action="store", dest="out_fasta",
                    help="output fasta file")
args = parser.parse_args()


#sequence_dict = SeqIO.index_db("temp_index.idx", [args.in_fasta], format="fasta")
sequence_dict = SeqIO.to_dict(SeqIO.parse(args.in_fasta, format="fasta"))
annotated_dict = {}
with open(args.in_gff, "r") as gff_fd:
    for record in GFF.parse(gff_fd, base_dict=sequence_dict):
        annotated_dict[record.id] = record

#print(annotated_dict['2R'].features[25])
with open(args.out_fasta, "w") as out_fd:
    for record in annotated_dict:
        for feature in annotated_dict[record].features:
            #print(feature.qualifiers)
            feature_location = "%s:%s-%s:%s" % (record, feature.location.start,
                                             feature.location.end, feature.location.strand)

            feature_id = ",".join(feature.qualifiers["Parent"]) if "Parent" in feature.qualifiers \
                           else ",".join(feature.qualifiers["ID"]) if "ID" in feature.qualifiers else "."
            feature_name = ",".join(feature.qualifiers["Name"]) if "Name" in feature.qualifiers else "."
            feature_seq = feature.extract(annotated_dict[record].seq)
            out_fd.write(">%s|%s|%s\n" % (feature_location, feature_id, feature_name))
            out_fd.write(str(feature_seq) + "\n")
#os.system("rm temp_index.idx")