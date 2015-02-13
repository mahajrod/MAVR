#!/usr/bin/env python
__author__ = 'mahajrod'

import os
import numpy as np

from Parsers.VCF import CollectionVCF
from Parsers.CCF import CollectionCCF
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.stats import chi2_contingency
from math import sqrt


def phi_coefficient_correlation(table):
    a, b, c, d = table[0][0], table[0][1], table[1][0], table[1][1]
    return (a*d - b*c) / sqrt((a+c)*(b+d)*(a+b)*(c+d))


def homogeneity_plot(sample_set_names_list, plot_file_prefix):
    rcParams.update({'font.size': 10})
    letter_list_part1 = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]
    feature_type_list = ["all", "pre_5'_UTR", "5'_UTR", "CDS", "3'_UTR"]

    plt.figure(1, dpi=300, figsize=(5*3, 3*len(sample_set_names_list)))
    index = 1
    for sample, letter in zip(sample_set_names_list, letter_list_part1):
        findex = 1
        for feature_type in feature_type_list:
            vcf_file = "%s_good.vcf" % sample if feature_type == "all" \
                else "%s_pre_UTR_variants_only_intergenic_l_300.vcf" % sample \
                if feature_type == "pre_5'_UTR" \
                else "%s_%s_variants.vcf" % (sample, feature_type)
            collection = CollectionVCF(from_file=True, vcf_file=vcf_file)
            sample_data = collection.count_strandness("%s_%s_variants" % (sample, feature_type))
            ax = plt.subplot(len(sample_set_names_list), len(feature_type_list), index)
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
            if index > len(feature_type_list) * (len(sample_set_names_list) - 1):
                plt.xlabel('Strandness')

            if findex == 1:
                plt.ylabel('N of SNV')
                plt.text(-0.4, 0.5, sample, rotation=90, fontweight="bold", transform=ax.transAxes, fontsize=14,
                         horizontalalignment='center',
                         verticalalignment='center')
            #plt.title("%s%i. %s (%i SNV)\np=%.3f, phi=%.3f" % (letter, findex, sample, len(collection), p_value, phi),
            #          fontweight='bold')
            #plt.title("%s%i. %s (%s)\np=%.3f, phi=%.3f" % (letter, findex, sample, feature_type, p_value, phi),
            #          fontweight='bold', fontsize=6)
            title = "%s%i.p=%.3f, phi=%.3f" % (letter, findex, p_value, phi) if p_value >= 0.001 \
                else "%s%i.p=%.2e, phi=%.3f" % (letter, findex, p_value, phi)

            plt.title(title, fontweight='bold', fontsize=10)
            plt.xticks(points + bar_width, ('None', '+', '-', 'Both'))
            if findex == len(feature_type_list):
                plt.legend()
            if index <= len(feature_type_list):
                plt.text(0.5, 1.15, feature_type, rotation=0, fontweight="bold", transform=ax.transAxes, fontsize=14,
                         horizontalalignment='center',
                         verticalalignment='center')
            #plt.suptitle("Strandness histograms", fontweight="bold", fontsize=20)
            findex += 1
            index += 1
    #plt.tight_layout()
    plt.savefig("%s.pdf" % plot_file_prefix)
    plt.savefig("%s.svg" % plot_file_prefix)
    plt.savefig("%s.eps" % plot_file_prefix)
    plt.close()

if __name__ == "__main__":
    workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/clusters/all/pre_UTR_strandness/"

    os.chdir(workdir)
    sample_set_names_list_6d = [#"PmCDA1_3d",
                                #"PmCDA1_sub1_3d",
                                "PmCDA1_6d",
                                "PmCDA1_sub1_6d"
                                ]

    sample_set_names_list_3d = [
                                "PmCDA1_3d",
                                "PmCDA1_sub1_3d",
                                #"PmCDA1_6d",
                                #"PmCDA1_sub1_6d"
                                ]

    homogeneity_plot(sample_set_names_list_3d, "PmCDA1_3d_strandness")
    homogeneity_plot(sample_set_names_list_6d, "PmCDA1_6d_strandness")


