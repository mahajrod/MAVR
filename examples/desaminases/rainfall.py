#!/usr/bin/env python

import os
from Parser.VCF import CollectionVCF, ReferenceGenome
from Parser.GFF import CollectionGFF
from BCBio import GFF


workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/filtered/"

sample_set_names_list   =   [#"PmCDA1_3d",
                             "HAP",
                             #"PmCDA1_sub1_3d",
                             #"PmCDA1_6d",
                             "HAP_sub1",
                             #"PmCDA1_sub1_6d",
                             #"A1_3d",
                             #"A1_6d",
                             #"A3G_3d",
                             #"AID_3d",
                             #"AID_6d"
                             ]
sample_suffix = "_adjusted_cluster_mutations.vcf"
rainfall_subdir = "rainfall_plot/"

reference = ReferenceGenome("/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/LAN210_v0.10m.fasta",
                                index_file="/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/LAN210_v0.10m.idx")


reference.find_gaps()

bad_regions_file = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/masked_regions/LAN210_v0.10m_masked_all_not_in_good_genes.gff"
bad_regions = CollectionGFF(input_file=bad_regions_file,
                                from_file=True)


bad_region_dict = {}
with open(bad_regions_file) as gff_fd:
    for record in GFF.parse(gff_fd):
        bad_region_dict[record.id] = record


os.chdir(workdir)
for sample_set in sample_set_names_list:
    mutations = CollectionVCF(from_file=True, vcf_file=sample_set + sample_suffix)
    mutations.rainfall_plot("%s_adjusted_cluster_mutations" % (sample_set), ref_genome=reference, draw_gaps=True,
                                        masked_regions=bad_region_dict, suptitle="%s" % sample_set, figsize=(36, 45))