#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os
from scipy import stats
from collections import OrderedDict

from CustomCollections.GeneralCollections import TwoLvlDict


def read_data(filename):
    data_dict = OrderedDict()
    with open(filename, "r") as in_fd:
        for line in in_fd:
            tmp = line.strip().split("\t")
            if tmp[0] == "YCL066W" or tmp[0] == "YCR097W_with_introns":
                data_dict[tmp[0]] = float(tmp[1])
    return data_dict


def calculate_chi_squared(data_dict, proportion_list):
    exp_value_list = [data_dict["YCL066W"], data_dict["YCR097W_with_introns"]]
    total = sum(exp_value_list)
    theor_value_list = []
    for proportion in proportion_list:
        theor_value_list.append(total * float(proportion) / sum(proportion_list))

    print exp_value_list
    print theor_value_list

    return stats.chisquare(exp_value_list, theor_value_list)


def get_results(samples_list, data_type):
    results = TwoLvlDict()

    for sample in samples_list:
        results[sample] = OrderedDict()
        filename = "%s/all_reads/%s_all_%s_coverage.tab" % (sample, sample, data_type)
        data = read_data(filename)
        if not data:
            print sample
            continue
        print sample
        for gene in data:
            results[sample][gene] = data[gene]

        for proportions, name in zip([[1, 2], [2, 1], [1, 1]], ["1:2", "2:1", "1:1"]):
            chi_results = calculate_chi_squared(data, proportions)
            print name
            results[sample][name + " Chi"] = chi_results[0]
            results[sample][name + " p-value"] = chi_results[1]
            print  chi_results
    return results


work_dir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/mating_type_detection"

os.chdir(work_dir)

temp = os.listdir("./")

samples_list = []
for filename in temp:
    if filename[0] == "N":
        samples_list.append(filename)

mean_results = get_results(samples_list, "mean")
median_results = get_results(samples_list, "median")

mean_results.write("all_mean_results.tab")
median_results.write("all_median_results.tab")