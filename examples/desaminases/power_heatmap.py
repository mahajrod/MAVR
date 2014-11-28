#!/usr/bin/env python

import os
import numpy as np

from Parser.VCF import ReferenceGenome, CollectionVCF, ref_alt_variants
from Parser.CCF import CollectionCCF
import matplotlib.pyplot as plt
import matplotlib as mpl

if __name__ == "__main__":
    workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/"

    sample_set_names_list = ["PmCDA1_3d",
                             "HAP",
                             "PmCDA1_sub1_3d",
                             "PmCDA1_6d",
                             "HAP_sub1",
                             "PmCDA1_sub1_6d",
    ]
    """
    sample_set_names_list = [#"HAP",
                             "HAP_sub1"]
    """

    clustering_dir = "clustering/"
    rainfall_dir = "rainfall"
    distance_threshold = 1000
    os.chdir(workdir)
    samples_list = sorted(os.listdir("."))
    min_cluster_size = 3
    # bad_regions_file = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/masked_regions/LAN210_v0.10m_masked_all.gff"


    annotations_dict = {}
    annotation_synonym_dict = {"three_prime_UTR": "3'_UTR",
                               "five_prime_UTR": "5'_UTR",
                               "snoRNA": "ncRNA",
                               "snRNA": "ncRNA"
    }
    clusters_dict = {}
    clusters_dict_3 = {}
    heatmap_dir = "heatmaps/"
    data = {}
    data_3 = {}
    y_name = "Power(Size/Median)"
    x_names_list = ["Size", "Median", "Homogeneity"]
    """
    for sample_set_name in sample_set_names_list:
        print("Handling %s" % sample_set_name)

        os.chdir(workdir)
        #os.system("mkdir -p %s" % sample_set_name)
        os.chdir(sample_set_name)
        #os.system("mkdir -p %s %s" % (clustering_dir, rainfall_dir))
        #os.system("pwd")
        clusters_dict[sample_set_name] = CollectionCCF(from_file=True, input_file="%s%s_adjusted_size_3+.ccf" %
                                                                                  (clustering_dir,
                                                                                   sample_set_name))

        data[sample_set_name] = clusters_dict[sample_set_name].get_data_for_stat(additional_data=("Median", "Power")) \
            if "HAP" in sample_set_name \
            else clusters_dict[sample_set_name].get_data_for_stat(additional_data=("Median", "Power", "Homogeneity"))

        data[sample_set_name] = [data[sample_set_name][:, i]
                                 for i in range(0, 4 if "HAP" in sample_set_name else 5)]


    os.chdir(workdir)
    os.system("mkdir -p %s" % heatmap_dir)
    os.chdir(heatmap_dir)

    x_names_dict = {"Length": 0, "Size": 1, "Median": 2, "Power": 3, "Homogeneity": 4}
    figure = plt.figure(1, dpi=150, figsize=(36, 18))
    subplot_dict = {}
    for i in range(0, len(x_names_list)):
        subplot_dict[x_names_list[i]] = []
        for j in range(0, len(sample_set_names_list)):
            if i == 2 and "HAP" in sample_set_names_list[j]:
                continue

            if i == 0 and j == 0:
                subplot_dict[x_names_list[i]].append(plt.subplot(3, 6, i * 6 + j + 1))
            else:
                if j == 0:
                    subplot_dict[x_names_list[i]].append(plt.subplot(3, 6, i * 6 + j + 1,
                                                                 sharey=subplot_dict[x_names_list[0]][0]))
                else:
                    subplot_dict[x_names_list[i]].append(plt.subplot(3, 6, i * 6 + j + 1,
                                                                 sharey=subplot_dict[x_names_list[0]][0],
                                                                 sharex=subplot_dict[x_names_list[i]][0]))
            x_data_index = x_names_dict[x_names_list[i]]
            x_data = data[sample_set_names_list[j]][x_data_index]
            y_data = data[sample_set_names_list[j]][3]

            max_y = max(y_data)
            max_x = max(x_data)
            n_y = int(max_y / 0.05) + 1

            if x_names_list[i] == "Size":
                n_x_bins = np.linspace(3, 30, 27 + 1)
                #print(n_x_bins)
                plt.xlim(xmax=30)
            elif x_names_list[i] == "Median":
                plt.xlim(xmax=1500)
                n_x_bins = np.linspace(0, 800, 16 + 1)
            elif x_names_list[i] == "Homogeneity":
                plt.xlim(xmax=1.0)
                n_x_bins = np.linspace(0.5, 1.0, 10 + 1)
            plt.xlabel(x_names_list[i])
            plt.ylabel("Power(Size/Median dist.)")
            plt.title(sample_set_names_list[j] + " (%i clusters)" % len(y_data))

            n_y_bins = np.linspace(0, 0.05 * n_y, n_y + 1)

            #print(i, j, sample_set_names_list[j])
            #plt.
            plt.hist2d(x_data, y_data, (n_x_bins, n_y_bins), cmin=1)
            plt.ylim(ymax=3)
            plt.colorbar()

    #figure.subplots_adjust(right=0.8)
    #cbar_ax = figure.add_axes([0.85, 0.15, 0.05, 0.7])
    #plt.colorbar(cax=cbar_ax)

    plt.savefig("heatmap_3+.svg")
    plt.close()
    """
    """
    for sample_set_name in sample_set_names_list:
        print("Handling %s" % sample_set_name)

        os.chdir(workdir)
        #os.system("mkdir -p %s" % sample_set_name)
        os.chdir(sample_set_name)
        #os.system("mkdir -p %s %s" % (clustering_dir, rainfall_dir))
        #os.system("pwd")
        clusters_dict[sample_set_name] = CollectionCCF(from_file=True, input_file="%s%s_adjusted_size_3+_power_0.05+.ccf" %
                                                                                  (clustering_dir,
                                                                                   sample_set_name))

        data[sample_set_name] = clusters_dict[sample_set_name].get_data_for_stat(additional_data=("Median", "Power")) \
            if "HAP" in sample_set_name \
            else clusters_dict[sample_set_name].get_data_for_stat(additional_data=("Median", "Power", "Homogeneity"))

        data[sample_set_name] = [data[sample_set_name][:, i]
                                 for i in range(0, 4 if "HAP" in sample_set_name else 5)]


    os.chdir(workdir)
    os.system("mkdir -p %s" % heatmap_dir)
    os.chdir(heatmap_dir)

    x_names_dict = {"Length": 0, "Size": 1, "Median": 2, "Power": 3, "Homogeneity": 4}
    figure = plt.figure(1, dpi=150, figsize=(36, 18))
    subplot_dict = {}
    for i in range(0, len(x_names_list)):
        subplot_dict[x_names_list[i]] = []
        for j in range(0, len(sample_set_names_list)):
            if i == 2 and "HAP" in sample_set_names_list[j]:
                continue

            if i == 0 and j == 0:
                subplot_dict[x_names_list[i]].append(plt.subplot(3, 6, i * 6 + j + 1))
            else:
                if j == 0:
                    subplot_dict[x_names_list[i]].append(plt.subplot(3, 6, i * 6 + j + 1,
                                                                 sharey=subplot_dict[x_names_list[0]][0]))
                else:
                    subplot_dict[x_names_list[i]].append(plt.subplot(3, 6, i * 6 + j + 1,
                                                                 sharey=subplot_dict[x_names_list[0]][0],
                                                                 sharex=subplot_dict[x_names_list[i]][0]))
            x_data_index = x_names_dict[x_names_list[i]]
            x_data = data[sample_set_names_list[j]][x_data_index]
            y_data = data[sample_set_names_list[j]][3]

            max_y = max(y_data)
            max_x = max(x_data)


            if x_names_list[i] == "Size":
                n_x_bins = np.linspace(3, 30, 27 + 1)
                #print(n_x_bins)
                plt.xlim(xmax=30)
            elif x_names_list[i] == "Median":
                plt.xlim(xmax=300)
                n_x_bins = np.linspace(0, 200, 10 + 1)
            elif x_names_list[i] == "Homogeneity":
                plt.xlim(xmax=1.0)
                n_x_bins = np.linspace(0.5, 1.0, 10 + 1)
            plt.xlabel(x_names_list[i])
            plt.ylabel("Power(Size/Median dist.)")
            plt.title(sample_set_names_list[j] + " (%i clusters)" % len(y_data))
            n_y = int(max_y / 0.05) + 1
            n_y_bins = np.linspace(0, 0.05 * n_y, n_y + 1)

            #print(i, j, sample_set_names_list[j])
            #plt.
            plt.hist2d(x_data, y_data, (n_x_bins, n_y_bins), cmin=1)
            plt.ylim(ymax=3)
            plt.colorbar()

    #figure.subplots_adjust(right=0.8)
    #cbar_ax = figure.add_axes([0.85, 0.15, 0.05, 0.7])
    #plt.colorbar(cax=cbar_ax)

    plt.savefig("heatmap_3+_power_0.05+.svg")
    plt.savefig("heatmap_3+_power_0.05+.eps")
    plt.close()
    """

    for sample_set_name in sample_set_names_list:
        print("Handling %s" % sample_set_name)

        os.chdir(workdir)
        #os.system("mkdir -p %s" % sample_set_name)
        os.chdir(sample_set_name)
        #os.system("mkdir -p %s %s" % (clustering_dir, rainfall_dir))
        #os.system("pwd")
        clusters_dict[sample_set_name] = CollectionCCF(from_file=True, input_file="%s%s_adjusted_size_3+_power_0.1+.ccf" %
                                                                                  (clustering_dir,
                                                                                   sample_set_name))

        data[sample_set_name] = clusters_dict[sample_set_name].get_data_for_stat(additional_data=("Median", "Power")) \
            if "HAP" in sample_set_name \
            else clusters_dict[sample_set_name].get_data_for_stat(additional_data=("Median", "Power", "Homogeneity"))

        data[sample_set_name] = [data[sample_set_name][:, i]
                                 for i in range(0, 4 if "HAP" in sample_set_name else 5)]


    os.chdir(workdir)
    os.system("mkdir -p %s" % heatmap_dir)
    os.chdir(heatmap_dir)



    x_names_dict = {"Length": 0, "Size": 1, "Median": 2, "Power": 3, "Homogeneity": 4}
    figure = plt.figure(1, dpi=150, figsize=(36, 18))
    subplot_dict = {}
    for i in range(0, len(x_names_list)):
        subplot_dict[x_names_list[i]] = []
        for j in range(0, len(sample_set_names_list)):
            if i == 2 and "HAP" in sample_set_names_list[j]:
                continue

            if i == 0 and j == 0:
                subplot_dict[x_names_list[i]].append(plt.subplot(3, 6, i * 6 + j + 1))
            else:
                if j == 0:
                    subplot_dict[x_names_list[i]].append(plt.subplot(3, 6, i * 6 + j + 1,
                                                                 sharey=subplot_dict[x_names_list[0]][0]))
                else:
                    subplot_dict[x_names_list[i]].append(plt.subplot(3, 6, i * 6 + j + 1,
                                                                 sharey=subplot_dict[x_names_list[0]][0],
                                                                 sharex=subplot_dict[x_names_list[i]][0]))
            x_data_index = x_names_dict[x_names_list[i]]
            x_data = data[sample_set_names_list[j]][x_data_index]
            y_data = data[sample_set_names_list[j]][3]

            max_y = max(y_data)
            max_x = max(x_data)

            if x_names_list[i] == "Size":
                n_x_bins = np.linspace(3, 30, 27 + 1)
                #print(n_x_bins)
                plt.xlim(xmax=30)
            elif x_names_list[i] == "Median":
                plt.xlim(xmax=300)
                n_x_bins = np.linspace(0, 200, 10 + 1)
            elif x_names_list[i] == "Homogeneity":
                plt.xlim(xmax=1.0)
                n_x_bins = np.linspace(0.5, 1.0, 10 + 1)
            plt.xlabel(x_names_list[i])
            plt.ylabel("Power(Size/Median dist.)")
            plt.title(sample_set_names_list[j] + " (%i clusters)" % len(y_data))
            n_y = int(max_y / 0.05) + 1
            n_y_bins = np.linspace(0, 0.05 * n_y, n_y + 1)

            #print(i, j, sample_set_names_list[j])
            #plt.
            plt.hist2d(x_data, y_data, (n_x_bins, n_y_bins), cmin=1)
            plt.ylim(ymax=3)
            plt.colorbar()

    #figure.subplots_adjust(right=0.8)
    #cbar_ax = figure.add_axes([0.85, 0.15, 0.05, 0.7])
    #plt.colorbar(cax=cbar_ax)

    plt.savefig("heatmap_3+_power_0.1+.svg")
    plt.savefig("heatmap_3+_power_0.1+.eps")
    plt.close()


