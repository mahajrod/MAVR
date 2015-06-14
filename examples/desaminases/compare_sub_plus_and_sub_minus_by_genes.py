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


def hypergeom(m, n, n1, n2, p=False):

    if m <= 0: return 1.0
    mmin = m - 1
    mmax = min(n1, n2)

    """
    return ss.hypergeom.cdf(mmax, n, n1, n2),\
           ss.hypergeom.cdf(mmin, n, n1, n2),\
           ss.hypergeom.cdf(mmax, n, n1, n2) - ss.hypergeom.cdf(mmin, n, n1, n2)

    return 1.0 - ss.hypergeom.sf(mmax, n, n1, n2),\
           1.0 - ss.hypergeom.sf(mmin, n, n1, n2),\
           -ss.hypergeom.sf(mmax, n, n1, n2) + ss.hypergeom.sf(mmin, n, n1, n2)
    """
    return -ss.hypergeom.sf(mmax, n, n1, n2) + ss.hypergeom.sf(mmin, n, n1, n2)


def get_intersection_length(start1, end1, start2, end2):
    if start1 - end2 > 0 or start2 - end1 > 0:
        return 0
    start_shift = start1 - start2
    start_coef_shift = 0 if start_shift < 0 else 1
    end_shift = end1 - end2
    end_coef_shift = 0 if end_shift > 0 else 1

    return (end2 - start2 + 1) - start_coef_shift * start_shift + end_coef_shift * end_shift

overlap_clusters_percent = TwoLvlDict({})

totaly_genes = 6074
test_fd = open("probability.t", "w")
test_fd.write("#size\tpower\ttotal\tPmCDA1_3d\tPmCDA1_sub_3d\tintersection\tp-value\n")
print([float(f) / float(100) for f in range(1, 11)])
for size in range(3, 11):
    overlap_clusters_percent[size] = {}
    for power in [float(f) / float(100) for f in range(1, 11)]:
        PmCDA1_3d_clusters = CollectionCCF(from_file=True, input_file="/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/clusters/%i/%.2f/PmCDA1_3d_size_%i+_power_%.2f+_good.ccf" % (size, power, size, power))
        PmCDA1_3d_sub_clusters = CollectionCCF(from_file=True, input_file="/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/clusters/%i/%.2f/PmCDA1_sub1_3d_size_%i+_power_%.2f+_good.ccf" % (size, power, size, power))

        #cluster_3d_dict = OrderedDict({})

        cluster_3d_set = set([])
        cluster_3d_sub_set = set([])

        for cluster_3d in PmCDA1_3d_clusters:
            for variant in cluster_3d:
                if "Genes" in variant.info_dict:
                    for gene in variant.info_dict["Genes"]:
                        cluster_3d_set.add(gene)

        for cluster_3d_sub in PmCDA1_3d_sub_clusters:
            for variant in cluster_3d_sub:
                if "Genes" in variant.info_dict:
                    for gene in variant.info_dict["Genes"]:
                        cluster_3d_sub_set.add(gene)


        print ("PmCDA1 3d : %i" % len(cluster_3d_set))
        print ("PmCDA1 3d sub: %i" % len(cluster_3d_sub_set))
        intersection = cluster_3d_set & cluster_3d_sub_set
        print("Intersection: % i" % len(intersection))

        n_intersection_genes = len(intersection)
        n_cluster_3d_genes = len(cluster_3d_set)
        n_cluster_3d_sub_genes = len(cluster_3d_sub_set)


        #print intersection
        #print cluster_3d_set
        #print cluster_3d_sub_set

        overlap_clusters_percent[size][power] = 100 * float(len(intersection))/float(len(cluster_3d_set))

        p_value = hypergeom(n_intersection_genes, totaly_genes, n_cluster_3d_genes, n_cluster_3d_sub_genes)
        test_fd.write("%i\t%.2f\t%i\t%i\t%i\t%i\t%e\n" % (size, power, totaly_genes, n_cluster_3d_genes, n_cluster_3d_sub_genes, n_intersection_genes, p_value))

overlap_clusters_percent.write("overlap_clusters_percent_genes.t")

test_fd.close()