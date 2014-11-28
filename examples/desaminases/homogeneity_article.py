#!/usr/bin/env python
__author__ = 'mahajrod'

import os
import numpy as np

from Parser.VCF import CollectionVCF
from Parser.CCF import CollectionCCF
import matplotlib.pyplot as plt
from matplotlib import rcParams

if __name__ == "__main__":
    workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/clusters/all/all/"
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
        collection = CollectionVCF(from_file=True, vcf_file=sample + "_good.vcf")
        sample_data = collection.count_strandness(sample + "_good_strandness")
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

        index += 1
    plt.savefig("PmCDA1_good_strandness.svg")
    plt.savefig("PmCDA1_good_strandness.eps")
    plt.close()

    """
    sample_set_names_list = ["PmCDA1_3d",
                             #"HAP",
                             "PmCDA1_sub1_3d",
                             "PmCDA1_6d",
                             #"HAP_sub1",
                             "PmCDA1_sub1_6d",
                             "A1_3d",
                             "A1_6d",
                             "A3G_3d",
                             "AID_3d",
                             "AID_6d"
                             ]
    for sample in sample_set_names_list:
        sample_data = CollectionVCF(from_file=True, vcf_file=sample + "_good.vcf").count_strandness(sample + "_good_strandness")
        plt.figure(1, dpi=300, figsize=(5, 5))
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
        plt.xlabel('Strandness')

        plt.ylabel('N of mutations')
        plt.title(sample)
        plt.xticks(points + bar_width, ('None', '+', '-', 'Both'))
        plt.legend()
        plt.savefig(sample + "_good_strandness.svg")
        plt.savefig(sample + "_good_strandness.eps")
        plt.close()
    """
    """
    sample_set_names_list = ["PmCDA1_3d",
                             "PmCDA1_sub1_3d",
                             "PmCDA1_6d",
                             "PmCDA1_sub1_6d",
                             "PmCDA1_3d",
                             "PmCDA1_sub1_3d",
                             "PmCDA1_6d",
                             "PmCDA1_sub1_6d"]

    for size in size_limits:
        os.chdir(workdir)
        for power in power_limits:
            os.chdir(workdir + size)
            os.chdir(power)
            data = {}
            plt.figure(1, dpi=300, figsize=(8, 4))
            index = 1
            subplot_list = []

            rcParams.update({'font.size': 5})

            for sample_set_name in sample_set_names_list:

                #print(index)
                prefix_filename = "%s_size_%s+_power_%s+_good" % (sample_set_name, size, power) if power != "all" \
                                  else "%s_size_%s+_good" % (sample_set_name, size)
                print("Handling %s" % prefix_filename)
                vcf_filename = prefix_filename + ".vcf"

                if index <= 4:
                    vcf_collection = CollectionVCF(from_file=True, vcf_file=vcf_filename)
                    mutation_count_file_prefix = "%s_size_%s+_power_%s+_good_strandness" % (sample_set_name, size, power)\
                                                 if power != "all" \
                                                 else "%s_size_%s+_good_strandness" % (sample_set_name, size)
                    mutation_data_dict = vcf_collection.count_strandness(mutation_count_file_prefix)
                    subplot_list.append(plt.subplot(2, len(sample_set_names_list)/2, index))

                    n_groups = 4
                    points = np.arange(n_groups)
                    bar_width = 0.35

                    C_values = mutation_data_dict["all"][0]
                    G_values = mutation_data_dict["all"][1]

                    rects1 = plt.bar(points, C_values, bar_width,
                                     color='b',
                                     label='C->T')

                    rects2 = plt.bar(points + bar_width, G_values, bar_width,
                                     color='g',
                                     label='G->A')
                    plt.xlabel('Strandness')
                    if index == 1:
                        plt.ylabel('N of mutations in clusters')
                    plt.title(sample_set_name)
                    plt.xticks(points + bar_width, ('None', '+', '-', 'Both'))
                    plt.legend()

                else:
                    ccf_filename = prefix_filename + ".ccf"
                    ccf_collection = CollectionCCF(from_file=True, input_file=ccf_filename)
                    cluster_data_dict = ccf_collection.get_data_for_stat(additional_data=["Homogeneity"])
                    if index == 5:
                        subplot_list.append(plt.subplot(2, len(sample_set_names_list)/2, index))
                    else:
                        subplot_list.append(plt.subplot(2, len(sample_set_names_list)/2, index, sharex=subplot_list[4]))

                    n_bins = np.linspace(0.5, 1.0, 4)

                    if index == 5:
                        plt.ylabel("N of clusters")
                    label = "Median: %f" % np.median(cluster_data_dict[:, 2])
                    plt.hist(cluster_data_dict[:, 2], bins=n_bins,
                             label=label)
                    plt.xlim(xmin=0.5, xmax=1.0)
                    plt.xlabel("Homogeneity\n(%s)" % label)
                    #plt.legend()
                index += 1
            plt.savefig("strandness_%s+_power_%s+.svg" % (size, power))
            plt.savefig("strandness_%s+_power_%s+.eps" % (size, power))
            plt.close()
    """
    """

    power_limits = ["%.2f" % (f / 100) for f in range(3, 11)] + ["all"]
    size_limits = [str(i) for i in range(5, 11)]

    clusters_dict = {}
    clusters_dict_3 = {}
    heatmap_dir = "heatmaps/"
    data = {}
    data_3 = {}
    y_name = "Power(Size/Median)"
    x_names_list = ["Size", "Homogeneity", "Power"]
    value_names_dict = {"Length": 0, "Size": 1, "Power": 2, "Homogeneity": 3}

    for size in size_limits:
        os.chdir(workdir)
        for power in power_limits:
            os.chdir(workdir + size)
            os.chdir(power)
            data = {}
            plt.figure(1, dpi=300, figsize=(8, 4))
            index = 1
            subplot_list = []

            rcParams.update({'font.size': 5})

            for sample_set_name in sample_set_names_list:

                #print(index)
                prefix_filename = "%s_size_%s+_power_%s+_good" % (sample_set_name, size, power) if power != "all" \
                                  else "%s_size_%s+_good" % (sample_set_name, size)
                print("Handling %s" % prefix_filename)
                vcf_filename = prefix_filename + ".vcf"

                if index <= 4:
                    vcf_collection = CollectionVCF(from_file=True, vcf_file=vcf_filename)
                    mutation_count_file_prefix = "%s_size_%s+_power_%s+_good_strandness" % (sample_set_name, size, power)\
                                                 if power != "all" \
                                                 else "%s_size_%s+_good_strandness" % (sample_set_name, size)
                    mutation_data_dict = vcf_collection.count_strandness(mutation_count_file_prefix)
                    subplot_list.append(plt.subplot(2, len(sample_set_names_list)/2, index))

                    n_groups = 4
                    points = np.arange(n_groups)
                    bar_width = 0.35

                    C_values = mutation_data_dict["all"][0]
                    G_values = mutation_data_dict["all"][1]

                    rects1 = plt.bar(points, C_values, bar_width,
                                     color='b',
                                     label='C->T')

                    rects2 = plt.bar(points + bar_width, G_values, bar_width,
                                     color='g',
                                     label='G->A')
                    plt.xlabel('Strandness')
                    if index == 1:
                        plt.ylabel('N of mutations in clusters')
                    plt.title(sample_set_name)
                    plt.xticks(points + bar_width, ('None', '+', '-', 'Both'))
                    plt.legend()

                else:
                    ccf_filename = prefix_filename + ".ccf"
                    ccf_collection = CollectionCCF(from_file=True, input_file=ccf_filename)
                    cluster_data_dict = ccf_collection.get_data_for_stat(additional_data=["Homogeneity"])
                    if index == 5:
                        subplot_list.append(plt.subplot(2, len(sample_set_names_list)/2, index))
                    else:
                        subplot_list.append(plt.subplot(2, len(sample_set_names_list)/2, index, sharex=subplot_list[4]))

                    n_bins = np.linspace(0.5, 1.0, 4)

                    if index == 5:
                        plt.ylabel("N of clusters")
                    label = "Median: %f" % np.median(cluster_data_dict[:, 2])
                    plt.hist(cluster_data_dict[:, 2], bins=n_bins,
                             label=label)
                    plt.xlim(xmin=0.5, xmax=1.0)
                    plt.xlabel("Homogeneity\n(%s)" % label)
                    #plt.legend()
                index += 1
            plt.savefig("strandness_%s+_power_%s+.svg" % (size, power))
            plt.savefig("strandness_%s+_power_%s+.eps" % (size, power))
            plt.close()
            """