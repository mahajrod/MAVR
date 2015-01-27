#!/usr/bin/env python2
__author__ = 'mahajrod'

import os
from collections import OrderedDict

from BCBio import GFF

from Bio import SeqIO, Seq
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from Parser.VCF import CollectionVCF


def count_mutations(mutations):
    count_dict = OrderedDict({})
    for record in mutations:
        try:
            if "Genes" in record.info_dict:
                for gene_name, gene_strand in zip(record.info_dict["Genes"], record.info_dict["Gstrand"]):
                    if gene_name not in count_dict:
                        count_dict[gene_name] = 0    # [0, 0, 0]
                    # mut_gene_strand - mutation strand with regard to gene strand
                    count_dict[gene_name] += 1
        except KeyError:
            print(record)
    return count_dict

sample_set_names_list = ["PmCDA1_3d",
                         "HAP",
                         "PmCDA1_sub1_3d",
                         "PmCDA1_6d",
                         "HAP_sub1",
                         "PmCDA1_sub1_6d",
                         ]
annotations_file = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/annotations/merged_annotations_Nagalakshmi_tranf_to_LAN210_v0.10m.gff3"
sequence_file = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/LAN210_v0.10m.fasta"
workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf"
sequence_dict = SeqIO.to_dict(SeqIO.parse(sequence_file, "fasta"))
#print(sequence_dict)

with open(annotations_file, "r") as in_fd:
    annotation_dict = dict([(record.id, record) for record in GFF.parse(in_fd)])

for sample_set in sample_set_names_list:
    os.chdir(workdir)
    os.chdir(sample_set)
    #muatations = CollectionVCF(from_file=True, vcf_file=)

tsv_fd = open("/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/promotor_regions.tsv", "w")
tsv_fd.write("gene_id\tsequence\n")
with open("/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/promotor_regions.fasta", "w") \
        as pr_fd:

    for record_id in annotation_dict:
        if record_id == "mt":
            continue
        for feature in annotation_dict[record_id].features:
            if feature.type != "gene":
                continue
            for subfeature in feature.sub_features:
                if subfeature.type == "five_prime_UTR":
                    break
            else:
                continue
            print(feature.strand)
            start = feature.location.start
            end = feature.location.end
            if feature.strand == +1:
                promotor_region = sequence_dict[record_id].seq[start-100:start]
            else:
                promotor_region = sequence_dict[record_id].seq[end:end+100].reverse_complement()

            print(feature.id, promotor_region, sequence_dict[record_id].seq[start:end])
            pr_fd.write(">%s\n%s\n" % (feature.id, str(promotor_region)))
            tsv_fd.write("%s\t%s\n" % (feature.id, str(promotor_region)))
tsv_fd.close()
#print("!!!!!!!!!!!!!!")
#print(sequence_dict)
    #record_dict = dict([(record.id, record) for record in GFF.parse(in_fd)])