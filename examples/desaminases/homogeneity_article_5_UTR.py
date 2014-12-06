#!/usr/bin/env python
__author__ = 'mahajrod'

import os
import numpy as np

from Parser.VCF import CollectionVCF
from Parser.CCF import CollectionCCF
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.stats import chi2_contingency
from math import sqrt


def phi_coefficient_correlation(table):
    a, b, c, d = table[0][0], table[0][1], table[1][0], table[1][1]
    return (a*d - b*c) / sqrt((a+c)*(b+d)*(a+b)*(c+d))

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

    feature_type_list = ["5'_UTR", "CDS", "3'_UTR"]

    for feature_type in feature_type_list:
        plt.figure(1, dpi=300, figsize=(6, 6))
        index = 1
        for sample, letter in zip(sample_set_names_list, letter_list_part1):
            collection = CollectionVCF(from_file=True, vcf_file="%s_%s_variants.vcf" % (sample, feature_type))
            sample_data = collection.count_strandness("%s_%s_variants" % (sample, feature_type))
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
            table = [C_values[1:-1], G_values[1:-1]]
            g, p_value, dof, expctd = chi2_contingency(table)
            phi = phi_coefficient_correlation(table)
            if index > 2:
                plt.xlabel('Strandness')
            if index == 1 or index == 3:
                plt.ylabel('N of SNV')
            plt.title("%s. %s (%i SNV)\np=%.3f, phi=%.3f" % (letter, sample, len(collection), p_value, phi),
                      fontweight='bold')
            plt.xticks(points + bar_width, ('None', '+', '-', 'Both'))
            plt.legend()
            plt.suptitle("%s_variants" % feature_type)
            index += 1
        plt.savefig("PmCDA1_%s_variants.svg" % (feature_type))
        plt.savefig("PmCDA1_%s_variants.eps" % (feature_type))
        plt.close()

