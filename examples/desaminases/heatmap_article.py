#!/usr/bin/env python

import os
import numpy as np

from Parsers.VCF import CollectionVCF
from Parsers.CCF import CollectionCCF
import matplotlib.pyplot as plt
from matplotlib import rcParams

if __name__ == "__main__":
    workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/clusters/"

    sample_set_names_list = ["HAP",
                             "PmCDA1_3d",
                             "PmCDA1_6d",
                             "HAP_sub1",
                             "PmCDA1_sub1_3d",
                             "PmCDA1_sub1_6d",
                             None,
                             "PmCDA1_3d",
                             "PmCDA1_6d",
                             None,
                             "PmCDA1_sub1_3d",
                             "PmCDA1_sub1_6d",

                             ]



    letter_list_part1 = ["A", "B", "C", "D", "E", "F", None, "G", "H", None, "I", "J"]
    os.chdir(workdir)

    power_limits = ["%.2f" % (f / 100) for f in range(3, 11)] + ["all"]
    size_limits = [str(i) for i in range(3, 11)]

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
            plt.figure(1, dpi=300, figsize=(6, 8))
            index = 1
            subplot_list = []

            rcParams.update({'font.size': 5})

            for sample_set_name, letter in zip(sample_set_names_list, letter_list_part1):
                #print(index)
                if index == 7 or index == 10:
                    index += 1
                    subplot_list.append(None)
                    continue
                file_name = "%s_size_%s+_power_%s+_good.ccf" % (sample_set_name, size, power) if power != "all" \
                            else "%s_size_%s+_good.ccf" % (sample_set_name, size)
                print("Handling %s" % file_name)
                data[sample_set_name] = {}
                if "HAP" in sample_set_name:
                    data_names = ["Size", "Power"]
                else:
                    data_names = ["Size", "Power", "Homogeneity"]

                tmp_data = CollectionCCF(from_file=True, input_file=file_name).get_data_for_stat(additional_data=data_names[1:])
                for data_name in data_names:
                    data[sample_set_name][data_name] = tmp_data[:, value_names_dict[data_name]]
                if index == 1:
                    subplot_list.append(plt.subplot(4, 3, index))
                elif index == 8:
                    subplot_list.append(plt.subplot(4, 3, index, sharex=subplot_list[0]))
                elif index >= 9:
                    subplot_list.append(plt.subplot(4, 3, index, sharex=subplot_list[7], sharey=subplot_list[7]))
                else:
                    subplot_list.append(plt.subplot(4, 3, index, sharex=subplot_list[0], sharey=subplot_list[0]))
                n_x_bins = np.linspace(int(size), 40, 40 - int(size) + 1)

                if index == 4 or index == 11 or index == 12:
                    plt.xlabel("Size")
                if 1 <= index <= 6:
                    n_y_bins = np.linspace(0, 6.0, int(6.0 / 0.05) + 1)
                else:
                    n_y_bins = np.linspace(0.5, 1.0, 11)
                if index == 1 or index == 4:
                    plt.ylabel("Power of cluster")
                elif index == 8 or index == 11:
                    plt.ylabel("Homogeneity")
                plt.title("%s. " % (letter) + sample_set_name + " (%i clusters)" % len(data[sample_set_name]["Size"]),
                          fontsize=5, fontweight='bold')
                counts, xedges, yedges, image = plt.hist2d(data[sample_set_name]["Size"],
                                                           data[sample_set_name]["Power"] if index <= 6 else data[sample_set_name]["Homogeneity"],
                                                           (n_x_bins, n_y_bins), cmin=1)

                plt.xlim(xmax=22)
                if index <= 6:
                    plt.ylim(ymax=2.5)
                else:
                    plt.ylim(ymax=1.0, ymin=0.5)
                #print(sample_set_name)
                #print(counts)

                max_counts = int(np.nanmax(counts))
                min_counts = 1
                print(index)
                #print(max_counts)
                #print(len(data[sample_set_name]["Size"]))

                #print(xedges)
                #print(yedges)
                cmap = plt.get_cmap('jet', max_counts)
                #cmap.set_under('gray')
                mappable = plt.cm.ScalarMappable(cmap=cmap)
                mappable.set_array([])
                mappable.set_clim(1, max_counts + 1)
                #mappable.set_array([])
                #mappable.set_clim(-0.5, ncolors+0.5)
                colorbar = plt.colorbar(mappable)
                #colorbar.set_ticks(np.linspace(1 + 0.5, max_counts + 0.5, max_counts), minor=True)
                decimal = int(max_counts / 10) + 1
                max_major_tick = max_counts - (max_counts % decimal)
                major_ticks = np.linspace(decimal + 0.5, max_major_tick + 0.5, int(max_major_tick / decimal))
                colorbar.set_ticks(major_ticks)
                colorbar.set_ticklabels(range(decimal, max_major_tick + decimal, decimal))

                if index <= 6:
                    max_size = max(data[sample_set_name]["Size"])
                    max_power = max(data[sample_set_name]["Power"])
                    x_bins = np.linspace(int(size), max_size + 1, max_size - int(size) + 2)
                    max_power_bin = (int(max_power / 0.01) + 1) * 0.01
                    min_power = float(power) if power != "all" else 0.0
                    y_bins = np.linspace(min_power, max_power_bin, int((max_power_bin - min_power)/0.01) + 1)
                    histogramm, xedges, yedges = np.histogram2d(data[sample_set_name]["Size"],
                                                                data[sample_set_name]["Power"],
                                                                bins=(x_bins, y_bins))

                    for array, array_name in zip([histogramm, xedges, yedges], ["histogramm", "xedges", "yedges"]):
                        np.savetxt("%s_%s_%s+_power_%s+.txt" % (array_name, sample_set_name, size, power),
                                      array, delimiter='\t', fmt="%s")

                index += 1

            plt.savefig("heatmap_%s+_power_%s+.svg" % (size, power))
            plt.savefig("heatmap_%s+_power_%s+.eps" % (size, power))
            plt.savefig("heatmap_%s+_power_%s+.pdf" % (size, power))
            plt.close()


