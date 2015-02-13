#!/usr/bin/env python2
__author__ = 'mahajrod'

import os
from Parsers.CCF import CollectionCCF
import matplotlib.pyplot as plt
import numpy as np

sample_set_names_list = ["PmCDA1_3d",
                             "HAP",
                             "PmCDA1_sub1_3d",
                             "PmCDA1_6d",
                             "HAP_sub1",
                             "PmCDA1_sub1_6d",
                             #"A1_3d",
                             #"A1_6d",
                             #"A3G_3d",
                             #"AID_3d",
                             #"AID_6d"
                             ]
power = "0.03"
size = "5"
bin_size = 50
workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/clusters/%s/%s/" % (size, power)
suffix = "_size_%s+_power_%s+_good.ccf" % (size, power)
hist_subfolder = "length_distribution/"
os.chdir(workdir)
os.system("mkdir -p %s" % hist_subfolder)

for extension in ".svg", ".eps", ".png", ".pdf":
    os.system("mkdir -p %s/%s" % (hist_subfolder, extension[1:]))
for sample_set in sample_set_names_list:
    collection = CollectionCCF(from_file=True, input_file=workdir+sample_set+suffix)
    length_data = collection.get_data_for_stat(additional_data=None)[:, 0]
    #print(stat_data)
    print(sample_set)
    #print(length_data)
    total = len(length_data)
    minimum = min(length_data)
    maximum = max(length_data)
    if len(length_data) == 0:
        continue
    plt.figure(1, figsize=(5, 5))
    plt.subplot(1, 1, 1)
    bins = np.linspace(1, 2500, 2500/bin_size+1)
    plt.hist(length_data, bins=bins, label="Min L: %i bp\nMax L: %i bp" % (minimum, maximum))
    plt.xlabel("Length of cluster")
    plt.ylabel("Number of clusters")
    plt.title("%s (%i clusters)" % (sample_set, total))
    plt.xlim(xmax=2500)
    plt.legend(fontsize=12)
    for extension in ".svg", ".eps", ".png", ".pdf":
        plt.savefig(workdir + hist_subfolder + extension[1:] + "/" + sample_set + "_length_distribution_size_%s+_power_%s+_binsize_%i" % (size, power, bin_size) + extension)
    plt.close()