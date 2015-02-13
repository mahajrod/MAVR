#!/usr/bin/env python
__author__ = 'mahajrod'

import os
from collections import OrderedDict
from Parsers.VCF import CollectionVCF
import numpy as np
import matplotlib.pyplot as plt

import pprint
from BCBio.GFF import GFFExaminer
from BCBio import GFF

if __name__ == "__main__":
    workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/clusters/all/all/"

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

    annotations = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/annotations/merged_annotations_Nagalakshmi_tranf_to_LAN210_v0.10m.gff3"
    with open(annotations, "r") as in_fd:
        record_dict = dict([(record.id, record) for record in GFF.parse(in_fd)])
    start_dict = OrderedDict({})
    end_dict = OrderedDict({})
    position_dict = OrderedDict({})
    length_dict = OrderedDict({})

    start_hist_dict = OrderedDict({})
    end_hist_dict = OrderedDict({})

    os.chdir(workdir)

    skip_genes_without_five_utr = False
    left = 300
    right = 300
    bin_width = 5
    bins = np.linspace(-left, right, ((left+right)/bin_width) + 1)
    normed = True
    max_start = 0
    max_end = 0
    for sample_set in sample_set_names_list:
        vcf_file = "%s_good.vcf" % sample_set
        #start_hist_prefix = "%s_start_hist_r_%i_l_%i" % (sample_set, right, left)
        #end_hist_prefix = "%s_end_hist_r_%i_l_%i" % (sample_set, right, left)
        #gene_variants = "%s_gene_variants_r_%i_l_%i.t" % (sample_set, right, left)
        variants = CollectionVCF(from_file=True, vcf_file=vcf_file)
        start_dict[sample_set], end_dict[sample_set], position_dict[sample_set] = \
            variants.variants_start_end(left, right, record_dict, skip_genes_without_five_utr=skip_genes_without_five_utr)
        length_dict[sample_set] = len(variants)
        #print(start_dict[sample_set])
        start_hist_dict[sample_set] = list(np.histogram(start_dict[sample_set], bins=bins))
        end_hist_dict[sample_set] = list(np.histogram(end_dict[sample_set], bins=bins))
        print(start_hist_dict[sample_set][0])
        if normed:
            start_hist_dict[sample_set][0] = start_hist_dict[sample_set][0].astype(np.float32, copy=False)
            end_hist_dict[sample_set][0] = end_hist_dict[sample_set][0].astype(np.float32, copy=False)
            start_hist_dict[sample_set][0] = start_hist_dict[sample_set][0] / length_dict[sample_set]
            end_hist_dict[sample_set][0] = end_hist_dict[sample_set][0] / length_dict[sample_set]
        print("Normed")
        print(start_hist_dict[sample_set][0])
        max_end = max(max_end, np.amax(end_hist_dict[sample_set][0]))
        max_start = max(max_start, np.amax(start_hist_dict[sample_set][0]))
    plt.figure(1, dpi=300, figsize=(16, 8*len(sample_set_names_list)))


    index = 1
    for sample_set in sample_set_names_list:
        plt.subplot(len(sample_set_names_list), 1, index)
        plt.bar(start_hist_dict[sample_set][1][:-1], start_hist_dict[sample_set][0], width=bin_width)
        plt.xlim(xmin=-left, xmax=right)
        #plt.ylim(ymax=max_start)
        plt.ylim(ymax=0.02)
        plt.axhline(0.02, color='y')
        plt.axhline(0.01, color='k')
        plt.axhline(0.005, color='r')
        plt.axhline(0.0025, color='g')
        plt.title(sample_set)
        index += 1
    skip = "_skipped_non_five_utr_genes" if skip_genes_without_five_utr else ""
    plt.savefig("all_DA_CDS_start_all_r_%i_l_%i_bin_width_%i%s.svg" % (right, left, bin_width, skip))
    plt.savefig("all_DA_CDS_start_all_r_%i_l_%i_bin_width_%i%s.eps" % (right, left, bin_width, skip))
    plt.close()