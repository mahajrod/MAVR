#!/usr/bin/env python2
import os

from BCBio import GFF

from Parser.VCF import ReferenceGenome, CollectionVCF, ref_alt_variants
from Parser.CCF import CollectionCCF
from Parser.GFF import CollectionGFF
if __name__ == "__main__":
    workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/"

    reference = ReferenceGenome("/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/LAN210_v0.10m.fasta",
                                index_file="/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/LAN210_v0.10m.idx")

    sample_set_names_list = ["PmCDA1_3d",
                             "PmCDA1_6d",
                             "A1_3d",
                             "A1_6d",
                             "A3G_3d",
                             "AID_3d",
                             "AID_6d",
                             "PmCDA1_sub1_3d",
                             "PmCDA1_sub1_6d",
                             "HAP",
                             "HAP_sub1"]
    """
    sample_set_names_list = [#"HAP",
                             "HAP_sub1"]
    """

    reference.find_gaps()
    os.chdir(workdir)
    samples_list = sorted(os.listdir("."))
    bad_regions_file = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/masked_regions/LAN210_v0.10m_masked_all.gff"
    bad_regions = CollectionGFF(input_file=bad_regions_file,
                                from_file=True)
    gff_file = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/annotations/merged_annotations_Nagalakshmi_tranf_to_LAN210_v0.10m.gff3"
    annotations_dict = {}
    annotation_synonym_dict = {"three_prime_UTR": "3'_UTR",
                               "five_prime_UTR": "5'_UTR",
                               "snoRNA": "ncRNA",
                               "snRNA": "ncRNA"
                               }
    annotation_black_list = ["region", "ARS", "long_terminal_repeat",
                                 "noncoding_exon", "intron", "repeat_region"]
    with open(gff_file) as gff_fd:
        for record in GFF.parse(gff_fd):
            annotations_dict[record.id] = record

    bad_region_dict = {}
    with open(bad_regions_file) as gff_fd:
        for record in GFF.parse(gff_fd):
            bad_region_dict[record.id] = record

    for sample_set_name in sample_set_names_list:
        print("Handling %s" % sample_set_name)
        os.chdir(workdir)
        mutations = CollectionVCF(vcf_file=sample_set_name + ".vcf",
                                  from_file=True)
        mutations.find_location(annotations_dict, use_synonym=True, synonym_dict=annotation_synonym_dict,
                                feature_type_black_list=annotation_black_list)
        mutations.write("%s_annotated.vcf" % sample_set_name)
