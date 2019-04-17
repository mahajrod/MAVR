#!/usr/bin/env python

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from Parsers.CCF import CollectionCCF

workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/clusters/"
all_files_subdir = "3/all/"
all_files_suffix = "_size_3+_good.ccf"
sample_set_names_list   =   ["PmCDA1_3d",
                             "HAP",
                             "PmCDA1_sub1_3d",
                             "PmCDA1_6d",
                             "HAP_sub1",
                             "PmCDA1_sub1_6d",

                             ]

samples_dir = workdir + all_files_subdir
power_limits = [f / 100 for f in range(1, 11)]
size_limits = [i for i in range(3, 11)]

os.chdir(workdir)
data_dict = {}
for sample in sample_set_names_list:
    data_dict[sample] = CollectionCCF(from_file=True, input_file=samples_dir + sample + all_files_suffix).get_data_for_stat(additional_data=["Power"])



figure = plt.figure(1, dpi=150, figsize=(18, 12))
for sample in sample_set_names_list:
    size_data = data_dict[sample][:, 1]
    max_size = max(size_data)
    bins_size_data = np.linspace(3, max_size, max_size + 1)
    power_data = data_dict[sample][:, 2]
    max_power = max(power_data)
    n_power_bins = int(max_power / 0.01) + 1
    bins_power_data = np.linspace(0, n_power_bins * 0.01, n_power_bins + 1)
    hist, xedges, yedges = np.histogram2d(size_data, power_data, bins=(bins_size_data, bins_power_data))
    ax = figure.add_subplot(3, 2, sample_set_names_list.index(sample) + 1, projection="3d")


    elements = (len(xedges) - 1) * (len(yedges) - 1)
    xpos, ypos = np.meshgrid(xedges[:-1]+0.25, yedges[:-1]+0.05)

    xpos = xpos.flatten()
    ypos = ypos.flatten()
    zpos = np.zeros(elements)
    dx = 0.5 * np.ones_like(zpos)
    dy = dx.copy()
    dz = hist.flatten()

    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color='b', zsort='average')

    #print (hist)
    #print (xedges)
    #print(yedges)


    plt.title(sample)
    plt.xlabel("Size")
    plt.ylabel("Power")
    plt.ylim(ymax=0.15)
plt.savefig("hist3d.svg")
plt.savefig("hist3d.eps")