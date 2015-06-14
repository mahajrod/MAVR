#!/usr/bin/env python

import os
import numpy as np
from collections import OrderedDict
from Parsers.VCF import ReferenceGenome, CollectionVCF, ref_alt_variants
from Parsers.CCF import CollectionCCF
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


if __name__ == "__main__":
    workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/"

    sample_set_names_list = ["PmCDA1_3d",
                             "PmCDA1_sub1_3d",
                             "PmCDA1_6d",
                             "PmCDA1_sub1_6d",
                             ]

    clustering_dir = "clustering/"
    os.chdir(workdir)

    annotation_synonym_dict = {"three_prime_UTR": "3'_UTR",
                               "five_prime_UTR": "5'_UTR",
                               "snoRNA": "ncRNA",
                               "snRNA": "ncRNA"
                               }

    suffix_list = ["_adjusted_size_3+.ccf",
                   "_adjusted_size_3+_power_0.05+.ccf",
                   "_adjusted_size_3+_power_0.1+.ccf"]

    homogeneity_dir = "homogeneity/"
    y_name = "N of clusters"

    os.system("mkdir -p %s" % homogeneity_dir)

    data = {}
    power_list = ["All", "Power >= 0.05", "Power >= 0.1"]
    for sample in sample_set_names_list:
        os.chdir(workdir + sample + "/" + clustering_dir)
        data[sample] = dict([(key, 0) for key in power_list])
        for suffix, name in zip(suffix_list, power_list):
            print(sample + suffix)
            clusters = CollectionCCF(from_file=True, input_file=sample + suffix)
            data[sample][name] = clusters.get_data_for_stat(additional_data=["Homogeneity", "Median", "Power"])

    parameters_dict = OrderedDict({"Length": 0,
                              "Size": 1,
                              "Homogeneity": 2,
                              "Median": 3,
                              "Power": 4})

    index = 1
    for parameter in parameters_dict:
        os.chdir(workdir + homogeneity_dir)
        plt.figure(index, dpi=150, figsize=(24, 18))
        for j in range(0, len(sample_set_names_list)):
            sample = sample_set_names_list[j]
            for i in range(0, len(power_list)):
                name = power_list[i]
                plt.subplot(3, 4, i * 4 + j + 1)
                plt.ylabel(y_name)
                plt.xlabel(parameter)
                tmp_data = data[sample][name][:, parameters_dict[parameter]]
                max_data = max(tmp_data)
                if parameter == "Size":
                    bins = np.linspace(3, max_data, max_data - 3 + 1)
                elif parameter == "Homogeneity":
                    bins = np.linspace(0.5, 1.0, 6)
                else:
                    bins = 20
                #print(parameter, bins)
                plt.hist(tmp_data, bins=bins, color="green")
                plt.title("%s (%s)" % (sample, name))
                if parameter == "Homogeneity":
                    plt.xlim(xmax=1)
                elif parameter == "Size":
                    plt.xlim(xmin=3, xmax=max_data)
                else:
                    plt.xlim(xmax=max_data)
        plt.savefig("%s_power_3x4.svg" % parameter)
        plt.savefig("%s_power_3x4.eps" % parameter)
        plt.close()
        index += 1













