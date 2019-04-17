#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from collections import OrderedDict
from RouToolPa.Routines import FileRoutines

def read_cluster_file(filename, with_counts=False):
    cluster_dict = OrderedDict()
    noname_family_index = 1
    with open(filename, "r") as in_fd:
        for line in in_fd:
            tmp = line.strip().split("\t")
            if with_counts:
                cluster_dict[tmp[0]] = (int(tmp[1]), set(tmp[2].split(",")))
            else:
                try:
                    tmp[1] = set(tmp[1].split(","))
                except IndexError:
                    tmp = ["Noname_fam_%i" % noname_family_index, set(tmp[0].split(","))]
                    noname_family_index += 1

                cluster_dict[tmp[0]] = (len(tmp[1]), tmp[1])

    return cluster_dict

parser = argparse.ArgumentParser()

parser.add_argument("-r", "--reference_file", action="store", dest="ref_file", required=True,
                    help="Input file with reference clusters.")
parser.add_argument("-c", "--file_to_check", action="store", dest="file_to_check", required=True,
                    help="File with clusters to compare with reference")
parser.add_argument("-u", "--ref_file_contains_gene_counts", action="store_true", dest="ref_with_counts",
                    default=False,
                    help="Reference file contains gene counts")
parser.add_argument("-n", "--check_file_contains_gene_counts", action="store_true", dest="check_with_counts",
                    default=False,
                    help="File to check contains gene counts")
parser.add_argument("-o", "--out_dir", action="store", dest="out_dir", default="compare_dir",
                    help="Output directory")

args = parser.parse_args()

FileRoutines.safe_mkdir(args.out_dir)

ref_clusters_dict = read_cluster_file(args.ref_file, with_counts=args.ref_with_counts)
check_clusters_dict = read_cluster_file(args.file_to_check, with_counts=args.check_with_counts)
totally_in_ref = len(ref_clusters_dict)
totally = len(check_clusters_dict)

synonym_file = "synonym.t"
contained_fully_in_file = "contained_fully_in.t"
contained_in_file = "contained_in.t"
include_file = "include.t"
all_file = "all.t"

synonym_dict = OrderedDict()
contained_fully_in_dict = OrderedDict()
contained_in_dict = OrderedDict()
include_dict = OrderedDict()

index = 1
for ref_cluster_id in ref_clusters_dict:
    for check_cluster_id in check_clusters_dict.keys():
        common_clusters = ref_clusters_dict[ref_cluster_id][1] & check_clusters_dict[check_cluster_id][1]
        common_clusters_len = len(common_clusters)
        if not common_clusters:
            continue
        if (common_clusters_len == ref_clusters_dict[ref_cluster_id][0]) and (common_clusters_len == check_clusters_dict[check_cluster_id][0]):
            synonym_dict[ref_cluster_id] = [check_cluster_id, ref_clusters_dict[ref_cluster_id][0], check_clusters_dict[check_cluster_id][0]]
            check_clusters_dict.pop(check_cluster_id, None)  # remove check_cluster_id from corresponding dict
            break
        if len(common_clusters) == ref_clusters_dict[ref_cluster_id][0]:
            # reference family is fully contained in checked cluster
            contained_fully_in_dict[ref_cluster_id] = [check_cluster_id, ref_clusters_dict[ref_cluster_id][0], check_clusters_dict[check_cluster_id][0]]
            break
        if len(common_clusters) == check_clusters_dict[check_cluster_id][0]:
            # reference family includes checked cluster
            if ref_cluster_id not in include_dict:
                include_dict[ref_cluster_id] = [check_cluster_id]
            else:
                include_dict[ref_cluster_id].append(check_cluster_id)
            continue

        # reference family is contained in checked clusters
        if ref_cluster_id not in contained_in_dict:
            contained_in_dict[ref_cluster_id] = [check_cluster_id]
        else:
            contained_in_dict[ref_cluster_id].append(check_cluster_id)
    else:

        if ref_cluster_id in include_dict:
            # checks in part of genes from reference cluster were not included in analysis
            if len(include_dict[ref_cluster_id]) == 1 and (ref_cluster_id not in contained_in_dict):
                #print ref_cluster_id
                synonym_dict[ref_cluster_id] = include_dict.pop(ref_cluster_id)
                #print synonym_dict[ref_cluster_id]
number_of_common_families = len(synonym_dict)
contained_fully_in_number = len(contained_fully_in_dict)
contained_in_number = len(contained_in_dict)
include_number = len(include_dict)

with open("%s/%s" % (args.out_dir, synonym_file), "w") as syn_fd:
    for fam_id in synonym_dict:
        #syn_fd.write("%s\t%s\t%i\t%i\n" % (fam_id, synonym_dict[fam_id][0], synonym_dict[fam_id][1], synonym_dict[fam_id][2]))
        syn_fd.write("%s\t%s\n" % (fam_id, synonym_dict[fam_id][0]))

with open("%s/%s" % (args.out_dir, all_file), "w") as syn_fd:
    for fam_id in ref_clusters_dict:
        #syn_fd.write("%s\t%s\t%i\t%i\n" % (fam_id, synonym_dict[fam_id][0], synonym_dict[fam_id][1], synonym_dict[fam_id][2]))
        #value = -1
        if fam_id in synonym_dict:
            value = synonym_dict[fam_id][0]
        elif fam_id in contained_fully_in_dict:
            # reference families fully contained in check families
            value = "C_%s" % contained_fully_in_dict[fam_id][0]
        elif fam_id in include_dict:
            # reference families that includes whole  several check families and in some cases parts of check families
            value = "I_%s" % ",".join(include_dict[fam_id])
            if fam_id in contained_in_dict:
                value += ";M_%s" % ",".join(contained_in_dict[fam_id])
        elif fam_id in contained_in_dict:
            value = "M_%s" % ",".join(contained_in_dict[fam_id])
        #if value == -1:
        #    value = "NF"
        syn_fd.write("%s\t%s\n" % (fam_id, value))


with open("%s/%s" % (args.out_dir, contained_fully_in_file), "w") as syn_fd:
    for fam_id in contained_fully_in_dict:
        #syn_fd.write("%s\t%s\t%i\t%i\n" % (fam_id, contained_fully_in_dict[fam_id][0], contained_fully_in_dict[fam_id][1], contained_fully_in_dict[fam_id][2]))
        syn_fd.write("%s\t%s\n" % (fam_id, contained_fully_in_dict[fam_id][0]))

with open("%s/%s" % (args.out_dir, contained_in_file), "w") as syn_fd:
    for fam_id in contained_in_dict:
        syn_fd.write("%s\t%s\n" % (fam_id, ",".join(contained_in_dict[fam_id])))

with open("%s/%s" % (args.out_dir, include_file), "w") as syn_fd:
    for fam_id in include_dict:
        syn_fd.write("%s\t%s\n" % (fam_id, ",".join(include_dict[fam_id])))

with open("%s/%s" % (args.out_dir, "stat.t"), "w") as syn_fd:
    syn_fd.write("Total_ref\t%i\n" % totally_in_ref)
    syn_fd.write("Total\t%i\n" % totally)
    syn_fd.write("Synonyms\t%i\nContains_fully_in\t%i\nContains_in\t%i\nIncludes_fully\t%i\n" % (number_of_common_families,
                                                                                    contained_fully_in_number,
                                                                                    contained_in_number,
                                                                                    include_number))