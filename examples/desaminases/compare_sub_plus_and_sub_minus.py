#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
from collections import OrderedDict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.stats as ss

from CustomCollections.GeneralCollections import TwoLvlDict
from Parsers.CCF import CollectionCCF


def get_intersection_length(start1, end1, start2, end2):
    if start1 - end2 > 0 or start2 - end1 > 0:
        return 0
    start_shift = start1 - start2
    start_coef_shift = 0 if start_shift < 0 else 1
    end_shift = end1 - end2
    end_coef_shift = 0 if end_shift > 0 else 1

    return (end2 - start2 + 1) - start_coef_shift * start_shift + end_coef_shift * end_shift

overlap_clusters_percent = TwoLvlDict({})
#size = 8
#power = 0.05
print([float(f) / float(100) for f in range(1, 11)])
for size in range(3, 11):
    overlap_clusters_percent[size] = {}
    for power in [float(f) / float(100) for f in range(1, 11)]:
        PmCDA1_3d_clusters = CollectionCCF(from_file=True, input_file="/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/clusters/%i/%.2f/PmCDA1_3d_size_%i+_power_%.2f+_good.ccf" % (size, power, size, power))

        PmCDA1_3d_sub_clusters = CollectionCCF(from_file=True, input_file="/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/clusters/%i/%.2f/PmCDA1_sub1_3d_size_%i+_power_%.2f+_good.ccf" % (size, power, size, power))
        PmCDA1_3d_clusters.write_gff("/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/clusters/%i/%.2f/PmCDA1_3d_size_%i+_power_%.2f+_good.gff" % (size, power, size, power))
        PmCDA1_3d_sub_clusters.write_gff("/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/clusters/%i/%.2f/PmCDA1_sub1_3d_size_%i+_power_%.2f+_good.gff" % (size, power, size, power))
        #cluster_3d_dict = OrderedDict({})

        cluster_3d_dict = TwoLvlDict({})

        for cluster_3d in PmCDA1_3d_clusters:
            cluster_3d_dict[cluster_3d.id] = OrderedDict({"length": cluster_3d.len,
                                                       "N of clusters": 0,
                                                       "length of clusters": [],
                                                       "intersection": [],
                                                       "intersection % of main cluster": [],
                                                       "interscection % of clusters": [],
                                                       "total_intersection": 0})
            for cluster_3d_sub in PmCDA1_3d_sub_clusters:
                if cluster_3d.chrom != cluster_3d_sub.chrom:
                    continue
                intersection = get_intersection_length(cluster_3d.start, cluster_3d.end, cluster_3d_sub.start, cluster_3d_sub.end)
                if intersection > 0:
                    cluster_3d_dict[cluster_3d.id]["N of clusters"] += 1
                    cluster_3d_dict[cluster_3d.id]["length of clusters"].append(cluster_3d_sub.len)
                    cluster_3d_dict[cluster_3d.id]["intersection"].append(intersection)
                    cluster_3d_dict[cluster_3d.id]["intersection % of main cluster"].append(intersection * 100/cluster_3d.len)
                    cluster_3d_dict[cluster_3d.id]["interscection % of clusters"].append(intersection * 100/cluster_3d_sub.len)

            cluster_3d_dict[cluster_3d.id]["total_intersection"] = sum(cluster_3d_dict[cluster_3d.id]["intersection % of main cluster"]) if cluster_3d_dict[cluster_3d.id]["intersection % of main cluster"] else 0

        cluster_3d_dict.write("intersection_PmCDA1_3d_sub_and_nonsub_%i+_%.2f+.t" % (size, power))

        total_intersection = [cluster_3d_dict[cluster_id]["total_intersection"] for cluster_id in cluster_3d_dict]
        print ("Total %i" % len(total_intersection))
        print("No intersection %i" % total_intersection.count(0))
        print("Intersection %i" % (len(total_intersection) - total_intersection.count(0)))
        figure = plt.figure(1, figsize=(5, 5), dpi=300)
        subplot = plt.subplot(1, 1, 1)
        plt.hist(total_intersection)
        plt.xlabel("% of intersection")
        plt.ylabel("N")
        plt.xlim(xmin=0, xmax=100)


        plt.savefig("intersection_PmCDA1_3d_sub_and_nonsub_%i+_%.2f+.svg" % (size, power))
        plt.close()
        overlap_clusters_percent[size][power] = 100 * float(len(total_intersection) - total_intersection.count(0))/float(len(total_intersection))

overlap_clusters_percent.write("overlap_clusters_percent.t")