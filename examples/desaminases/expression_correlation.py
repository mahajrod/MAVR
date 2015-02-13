#!/usr/bin/env python2

__author__ = 'mahajrod'

import os
from collections import OrderedDict
from copy import deepcopy

import numpy as np

from Bio import SeqIO
from BCBio import GFF

from Parsers.VCF import CollectionVCF
from Parsers.Cufflinks import CollectionFPKMTracking
from Parsers.GFF import CollectionGFF


def prepare_data(mutations, prefix, expression_data_dir, sequence_dict, annotation_dict):
    count_dict = OrderedDict({})
    for record in mutations:
        mutation_strand = mutation_strand_dict[record.ref]
        try:
            if "Genes" in record.info_dict:
                for gene_name, gene_strand in zip(record.info_dict["Genes"], record.info_dict["Gstrand"]):
                    if gene_name not in count_dict:
                        count_dict[gene_name] = np.zeros(4)     # [0, 0, 0]
                    # mut_gene_strand - mutation strand with regard to gene strand
                    mut_gene_strand = "P" if mutation_strand == gene_strand else "M"
                    count_dict[gene_name][values_names[mut_gene_strand]] += 1
        except KeyError:
            print(record)



    print("Totaly %s genes with mutations" % len(count_dict))
    os.system("mkdir -p %s" % expression_cor_dir)

    for directory in expression_samples_dir_list:
        os.chdir("%s%s/%s" % (workdir, sample_set_name, expression_cor_dir))
        os.system("mkdir -p %s" % directory)
        os.chdir(directory)
        gene_expresion_data = CollectionFPKMTracking(from_file=True, input_file="%s%sfiltered_genes.fpkm_tracking"
                                                                                % (expression_data_dir, directory))
        temp_dict = deepcopy(count_dict)

        for record in gene_expresion_data:
            gene_name = record.gene_id if record.gene_id != "" else \
                record.gene_short_name if record.gene_short_name != "" else None
            if gene_name is None:
                continue

            if gene_name in temp_dict:
                temp_dict[gene_name][values_names["Exp"]] = record.FPKM
            else:
                temp_dict[gene_name] = np.zeros(4)
                temp_dict[gene_name][values_names["Exp"]] = record.FPKM

            temp_dict[gene_name][values_names["Len"]] = record.length

        with open("%s_mut_expression.t" % prefix, "w") as out_fd:
            out_fd.write("gene\tLen\tP\tN\tExp\n")
            for gene_name in temp_dict:
                out_fd.write("%s\t%s\n" % (gene_name, "\t".join([str(x) for x in temp_dict[gene_name]])))


if __name__ == "__main__":
    workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/"

    sample_set_names_list = ["PmCDA1_3d",
                             "PmCDA1_6d",
                             "A1_3d",
                             "A1_6d",
                             "A3G_3d",
                             "AID_3d",
                             "AID_6d",
                             "PmCDA1_sub1_3d",
                             "PmCDA1_sub1_6d",
                             #"HAP",
                             #"HAP_sub1"
                             ]
    """
    expression_data_dir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/expression/Yeast_RNA_Seq/"
    expression_samples_dir_list = ["S1_Nagal_clout/",
                                   #"S1_Yassour_clout/",
                                   "S2_Nagal_clout/",
                                   #"S2_Yassour_clout/"
                                   ]
    expression_cor_dir = "expression_correlation/"
    """
    expression_data_dir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/expression/stranded_expression/"
    expression_samples_dir_list = ["S1_Nagal_clout/",
                                   "S1_Yassour_clout/",
                                   "S2_Nagal_clout/",
                                   "S2_Yassour_clout/",
                                   "S3_Nagal_clout/",
                                   "S3_Yassour_clout/",
                                   "S4_Nagal_clout/",
                                   "S4_Yassour_clout/",
                                   "S5_Nagal_clout/",
                                   "S5_Yassour_clout/",
                                   "S6_Nagal_clout/",
                                   "S6_Yassour_clout/"]
    expression_cor_dir = "stranded_expression_correlation/"


    annotations_file = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/annotations/merged_annotations_Nagalakshmi_tranf_to_LAN210_v0.10m.gff3"
    sequence_file = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/LAN210_v0.10m.fasta"
    workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf"
    sequence_dict = SeqIO.to_dict(SeqIO.parse(sequence_file, "fasta"))
    #print(sequence_dict)

    with open(annotations_file, "r") as in_fd:
        annotation_dict = dict([(record.id, record) for record in GFF.parse(in_fd)])



    mutation_strand_dict = {"C": "P", "G": "M"}

    #values_names = {"Len": 0, "P": 1, "P_dens": 2, "M": 3, "M_dens": 4, "Exp": 5}
    values_names = {"Len": 0, "P": 1, "M": 2, "Exp": 3}
    for sample_set_name in sample_set_names_list:
        print("Handling %s" % sample_set_name)
        os.chdir(workdir)

        os.chdir(sample_set_name)
        mutations = CollectionVCF(vcf_file="./clustering/%s_adjusted_cluster_mutations.vcf" % sample_set_name,
                                  from_file=True)
        mutations_large_clusters = CollectionVCF(vcf_file="./clustering/%s_adjusted_3+_cluster_mutations.vcf"
                                                          % sample_set_name,
                                                 from_file=True)
        prepare_data(mutations, "all_adjusted_cluster_mutations", expression_data_dir)
        prepare_data(mutations_large_clusters, "adjusted_3+_cluster_mutations", expression_data_dir)