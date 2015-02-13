#!/usr/bin/env python2
import os
import numpy as np
from collections import OrderedDict
from math import sqrt
import matplotlib.pyplot as plt

from Parsers.VCF import ReferenceGenome, CollectionVCF


def draw_bar_diagram(count_dict, out_file, figsize=(10, 10)):
    num_of_regions = len(count_dict)
    side = int(sqrt(num_of_regions))
    if side*side != num_of_regions:
        side += 1

    fig = plt.figure(1, dpi=200, figsize=figsize, facecolor="#D6D6D6")
    fig.suptitle("Strandness diagram", fontsize=16)
    sub_plot_dict = OrderedDict({})
    index = 1
    n_groups = 4
    points = np.arange(n_groups)
    bar_width = 0.35
    for chromosome in count_dict:
        C_values = count_dict[chromosome][0]
        G_values = count_dict[chromosome][1]
        sub_plot_dict[chromosome] = plt.subplot(side, side, index, axisbg="#D6D6D6")

        index += 1

        rects1 = plt.bar(points, C_values, bar_width,
                         color='b',
                         label='C')

        rects2 = plt.bar(points + bar_width, G_values, bar_width,
                         color='g',
                         label='G')
        plt.xlabel('Strandness')
        plt.ylabel('Counts')
        plt.title(chromosome)
        plt.xticks(points + bar_width, ('N', 'P', 'M', 'B'))
        plt.legend()

    plt.savefig(out_file)
    plt.close()


def count_strandness(mutations, prefix):
    count_dict = OrderedDict({})

    hor_coord_dict = {"C": 0, "G": 1}
    ver_coord_dict = {"N": 0, "P": 1, "M": 2, "B": 3}
    for record in mutations:
        if record.chrom not in count_dict:
            count_dict[record.chrom] = np.zeros((2, 4), dtype=int)
        count_dict[record.chrom][hor_coord_dict[record.ref]][ver_coord_dict[record.info_dict["Fstrand"][0]]] += 1

    draw_bar_diagram(count_dict, "%s_chromosome_strandnes_diagram.svg" % prefix, figsize=(20, 20))

    count_dict["all"] = sum(count_dict.values())

    draw_bar_diagram({"all": count_dict["all"]}, "%s_all_chromosome_strandnes_diagram.svg" % prefix, figsize=(5, 5))

    for chromosome in count_dict:
        with open("%s_%s.t" % (prefix, chromosome), "w") as out_fd:
            out_list = count_dict[chromosome].tolist()
            for index, name in zip(range(0, len(out_list)), ["C", "G"]):
                out_list[index].insert(0, name)
            out_list.insert(0, [".", "N", "P", "M", "B"])
            for string_list in out_list:
                out_fd.write("\t".join([str(x) for x in string_list]) + "\n")

    return(count_dict)


if __name__ == "__main__":
    workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/"

    reference = ReferenceGenome("/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/LAN210_v0.10m.fasta",
                                index_file="/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/LAN210_v0.10m.idx")

    sample_set_names_list = ["PmCDA1_3d",
                             "PmCDA1_6d",
                             "A1_3d",
                             "A1_6d",
                             "A3G_3d",
                             "AID_3d",
                             "AID_6d",
                             "PmCDA1_sub1_3d",
                             "PmCDA1_sub1_6d",
                             #"HAP",
                             #"HAP_sub1"
                             ]
    """
    sample_set_names_list = [#"HAP",
                             "HAP_sub1"]
    """

    reference.find_gaps()
    os.chdir(workdir)
    samples_list = sorted(os.listdir("."))

    for sample_set_name in sample_set_names_list:
        print("Handling %s" % sample_set_name)
        os.chdir(workdir)
        os.chdir(sample_set_name)
        mutations = CollectionVCF(vcf_file="./clustering/%s_adjusted_cluster_mutations.vcf" % sample_set_name,
                                  from_file=True)
        mutations_large_clusters = CollectionVCF(vcf_file="./clustering/%s_adjusted_3+_cluster_mutations.vcf"
                                                          % sample_set_name,
                                                 from_file=True)
        os.system("mkdir -p strandness")
        os.chdir("strandness")

        mut_count = count_strandness(mutations, "all_good_mut_strandness")
        mut_large_count = count_strandness(mutations_large_clusters, "3+_clusters_good_mut_strandness")


