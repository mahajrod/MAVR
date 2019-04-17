#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
from RouToolPa.Tools.GATK import FastaAlternateReferenceMaker
from RouToolPa.Routines.File import read_ids




ref_dir = "/home/mahajrod/Genetics/Projects/desaminases/data/LAN210_v0.10m/"
reference = ref_dir + "LAN210_v0.10m.fasta"
ref_annotations = ref_dir + "annotations/merged_annotations_Nagalakshmi_tranf_to_LAN210_v0.10m.gff3"
gene_ids_file = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/polymorphisms/gene.ids"
gatk_dir = "/home/mahajrod/Repositories/genetic/NGS_tools/GenomeAnalysisTK-3.2-0/"

work_dir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/"

gene_ids = read_ids(gene_ids_file)

sample_set_names_list = ["PmCDA1_3d",
                          "HAP",
                             "PmCDA1_sub1_3d",
                             "PmCDA1_6d",
                             "HAP_sub1",
                             "PmCDA1_sub1_6d",
                             "A1_3d",
                             "A1_6d",
                             "A3G_3d",
                             "AID_3d",
                             "AID_6d"
                             ]

for sample in sample_set_names_list:
    os.chdir(work_dir)
    os.chdir(sample)
    os.system("mkdir -p genome")
    variants_file_list = os.listdir("per_sample_vcf/")
    for vcf_file in variants_file_list:
        FastaAlternateReferenceMaker.correct_reference(gatk_dir, reference, "genome/%s.fasta" % vcf_file[:-4],
                                                       "per_sample_vcf/%s" % vcf_file)