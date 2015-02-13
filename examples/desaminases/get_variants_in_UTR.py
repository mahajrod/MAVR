#!/usr/bin/env python
__author__ = 'mahajrod'

import os
from Parsers.VCF import CollectionVCF


def variants_in_UTR(collection, feature_type="5'_UTR"):

    UTR_record_list = []
    for variant in collection:
        if feature_type in variant.info_dict["Ftype"]:
            UTR_record_list.append(variant)
    return CollectionVCF(from_file=False, record_list=UTR_record_list, metadata=collection.metadata, header=collection.header)

if __name__ == "__main__":
    workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/clusters/all/pre_UTR_strandness/"

    sample_set_names_list = ["PmCDA1_3d",
                             #"HAP",
                             "PmCDA1_sub1_3d",
                             "PmCDA1_6d",
                             #"HAP_sub1",
                             "PmCDA1_sub1_6d",
                             #"A1_3d",
                             #"A1_6d",
                             #"A3G_3d",
                             #"AID_3d",
                             #"AID_6d"
                             ]
    os.chdir(workdir)
    feature_type_list = ["5'_UTR", "CDS", "3'_UTR"]

    for sample_set in sample_set_names_list:
        vcf_file = "%s_good.vcf" % sample_set
        variants = CollectionVCF(from_file=True, vcf_file=vcf_file)
        for feature_type in feature_type_list:
            pre_UTR_variants = variants_in_UTR(variants, feature_type=feature_type)
            pre_UTR_variants.write("%s_%s_variants.vcf" % (sample_set, feature_type))
