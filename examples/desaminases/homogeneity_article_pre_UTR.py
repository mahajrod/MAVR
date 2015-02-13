#!/usr/bin/env python
__author__ = 'mahajrod'

import os
import numpy as np

from Parsers.VCF import CollectionVCF
from Parsers.CCF import CollectionCCF
import matplotlib.pyplot as plt
from matplotlib import rcParams

if __name__ == "__main__":
    workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/clusters/all/pre_UTR_strandness/"
    letter_list_part1 = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]
    os.chdir(workdir)
    sample_set_names_list = ["PmCDA1_3d",
                             "PmCDA1_sub1_3d",
                             "PmCDA1_6d",
                             "PmCDA1_sub1_6d"
                             ]
    rcParams.update({'font.size': 7})
    plt.figure(1, dpi=300, figsize=(6, 6))
    index = 1
    for sample, letter in zip(sample_set_names_list, letter_list_part1):
        collection = CollectionVCF(from_file=True, vcf_file=sample + "_pre_UTR_variants_only_intergenic_l_300.vcf")
        sample_data = collection.count_strandness(sample + "_pre_UTR_variants_only_intergenic_l_300_strandness")
        plt.subplot(2, 2, index)
        n_groups = 4
        points = np.arange(n_groups)
        bar_width = 0.35

        C_values = sample_data["all"][0]
        G_values = sample_data["all"][1]

        rects1 = plt.bar(points, C_values, bar_width,
                                     color='b',
                                     label='C->T')

        rects2 = plt.bar(points + bar_width, G_values, bar_width,
                                     color='g',
                                     label='G->A')
        if index > 2:
            plt.xlabel('Strandness')
        if index == 1 or index == 3:
            plt.ylabel('N of SNV')
        plt.title("%s. %s (%i SNV) " % (letter, sample, len(collection)), fontweight='bold')
        plt.xticks(points + bar_width, ('None', '+', '-', 'Both'))
        plt.legend()
        plt.suptitle("Pre_UTR_variants_only_intergenic_l_300")
        index += 1
    plt.savefig("PmCDA1_pre_UTR_variants_only_intergenic_l_300_strandness.svg")
    plt.savefig("PmCDA1_pre_UTR_variants_only_intergenic_l_300_strandness.eps")
    plt.close()

