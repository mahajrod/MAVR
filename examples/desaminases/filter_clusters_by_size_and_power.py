#!/usr/bin/env python

import os

from Parser.CCF import CollectionCCF
from Parser.VCF import CollectionVCF

from Collections.GeneralCollections import TwoLvlDict
from collections import OrderedDict

workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/clusters/"
all_files_subdir = "all/all/"
all_files_suffix = "_good.ccf"
sample_set_names_list   =   ["PmCDA1_3d",
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
power_limits = [f / 100 for f in range(1, 11)]
size_limits = [i for i in range(3, 11)]

os.chdir(workdir)
for sample_set in sample_set_names_list:
    stat_dict = TwoLvlDict(OrderedDict({}))
    print("Handling %s" % sample_set)
    all_clusters = CollectionCCF(from_file=True, input_file=workdir + all_files_subdir + sample_set + all_files_suffix)
    if "HAP" not in sample_set:
        all_clusters.check_strandness()
    for min_size in size_limits:
        stat_dict[min_size] = OrderedDict({})
        os.system("mkdir -p %i %i/all " % (min_size, min_size))
        above_size_clusters, below_size_clusters = all_clusters.filter_by_expression("record.size >= %i" % min_size)
        above_size_clusters.write("%i/all/%s_size_%i+%s" % (min_size, sample_set, min_size, all_files_suffix))
        stat_dict[min_size][0.00] = len(above_size_clusters)
        for min_power in power_limits:

            os.system("mkdir -p %i/%.2f" % (min_size, min_power))
            above_power_clusters, below_power_clusters = above_size_clusters.filter_by_expression("record.description['Power'] >= %f" % min_power)
            above_power_clusters.write("%i/%.2f/%s_size_%i+_power_%.2f+%s" % (min_size, min_power, sample_set, min_size, min_power, all_files_suffix))
            stat_dict[min_size][min_power] = len(above_power_clusters)

    stat_dict.write("%s_statistics.t" % sample_set)


