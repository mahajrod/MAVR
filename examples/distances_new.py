#!/usr/bin/env python
__author__ = 'mahajrod'


import os, sys
import matplotlib.pyplot as plt
import numpy as np
sys.path.append("../")

from Parse.ParseR import *
from SeqAnalysis.SeqAnalysis import get_taxonomy_from_genbank_files
from Pipelines.Phylogeny import get_distances


def ismember(a1, a2):
    """ Test whether items from a2 are in a1.

    This does the same thing as np.setmember1d, but works on
    non-unique arrays.

    Only a few (2-4) times slower than np.setmember1d, and a lot
    faster than [i in a2 for i in a1].

    An example that np.setmember1d gets wrong:

    >>> a1 = np.array([5,4,5,3,4,4,3,4,3,5,2,1,5,5])
    >>> a2 = [2,3,4]
    >>> mask = ismember(a1,a2)
    >>> a1[mask]
    array([4, 3, 4, 4, 3, 4, 3, 2])
    """
    a2 = set(a2)
    a1 = np.asarray(a1)
    ind = a1.argsort()
    a1 = a1[ind]
    mask  = []
    # need this bit because prev is not defined for first item
    item  = a1[0]
    if item in a2:
        mask.append(True)
        a2.remove(item)
    else:
        mask.append(False)
    prev = item
    # main loop
    for item in a1[1:]:
        if item == prev:
            mask.append(mask[-1])
        elif item in a2:
            mask.append(True)
            prev = item
            a2.remove(item)
        else:
            mask.append(False)
            prev = item
    # restore mask to original ordering of a1 and return
    mask = np.array(mask)
    return mask[ind.argsort()]

taxa_file = "/home/mahajrod/genetics/MH_selection/data/complete_genomes_13.05.2014/splited/random_genomes/Actinopterygii_Dipnoi/Actinopterygii_Dipnoi_random_genomes.gb"
genus_dict, family_dict, order_dict = get_taxonomy_from_genbank_files([taxa_file], "new_data.idx")

tree_dist_name = "/home/mahajrod/genetics/MH_selection/data/complete_genomes_13.05.2014/splited/random_genomes/Actinopterygii_Dipnoi/iteration_2/partition_finder/tree_dist.csv"
column_names, row_names, table = parse_R_table_file(tree_dist_name, separator=",")

genus_distances = []
family_distances = []
order_distances = []

genus_distances = get_distances(genus_dict, row_names, table)
order_distances = get_distances(order_dict, row_names, table)
family_distances = get_distances(family_dict, row_names, table)

np.savetxt("genus_distances.t", genus_distances)
np.savetxt("family_distances.t", family_distances)
np.savetxt("order_distances.t", order_distances)

plt.figure(1, dpi=400, figsize=(16, 8))
bins = np.linspace(0, 3.0, 100)

plt.hist([genus_distances, family_distances, order_distances],
         bins,
         alpha=0.5,
         label=['genus', 'family', "order"],
         color=["blue", "green", "red"])

plt.title("Pairwise evolutionary distances between species")
plt.legend(loc='upper right')
plt.xlabel("distances")
plt.ylabel("N of pairs")
plt.savefig("all.svg")
plt.figure(2, dpi=400, figsize=(16, 8))
bins = np.linspace(0, 3.0, 100)

order_indexes = np.nonzero(ismember(order_distances, family_distances))[0]
family_indexes = np.nonzero(ismember(family_distances, genus_distances))[0]

only_order_distances = np.delete(order_distances, order_indexes)
only_family_distances = np.delete(family_distances, family_indexes)

np.savetxt("only_family_distances.t", only_family_distances)
np.savetxt("only_order_distances.t", only_order_distances)

plt.hist([genus_distances, only_family_distances, only_order_distances],
         bins,
         alpha=0.5,
         label=['genus', 'family', "order"],
         color=["blue", "green", "red"])
plt.title("Pairwise evolutionary distances between species")
plt.legend(loc='upper right')
plt.xlabel("distances")
plt.ylabel("N of pairs")
plt.savefig("all_no_intersect.svg")


random_order_distances = np.random.choice(only_order_distances, len(genus_distances), replace=False)
random_family_distances = np.random.choice(only_family_distances, len(genus_distances), replace=False)

np.savetxt("random_only_family_distances.t", random_family_distances)
np.savetxt("random_only_order_distances.t", random_order_distances)

plt.figure(3, dpi=400, figsize=(12, 6))
bins = np.linspace(0, 3.0, 50)
plt.hist([genus_distances, random_family_distances, random_order_distances],
         bins,
         alpha=0.5,
         label=['genus', 'family', "order"],
         color=["blue", "green", "red"])
plt.title("Pairwise evolutionary distances between species")
plt.legend(loc='upper right')
plt.xlabel("distances")
plt.ylabel("N of pairs")
plt.savefig("all_no_intersect_random.svg")
plt.savefig("all_no_intersect_random.png")