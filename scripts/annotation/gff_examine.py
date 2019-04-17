#!/usr/bin/env python
__author__ = 'mahajrod'
import os
import argparse
import pprint
from collections import OrderedDict
import matplotlib



matplotlib.use('Agg')
os.environ['MPLCONFIGDIR'] = '/tmp/'
import matplotlib.pyplot as plt
plt.ioff()
import numpy as np

from BCBio.GFF import GFFExaminer
from BCBio import GFF
from RouToolPa.Collections.General import TwoLvlDict
from RouToolPa.Routines.Sequence import get_feature_lengths, get_total_feature_lengths, feature_lengths_collapse_records

parser = argparse.ArgumentParser()

parser.add_argument("-f", "--gff", action="store", dest="gff",
                    help="gff to examine")
parser.add_argument("-o", "--len_file", action="store", dest="len_file",
                    help="Output file for feature lengths", default=None)
parser.add_argument("-p", "--prefix", action="store",
                    dest="prefix",
                    help="Prefix of output files",
                    default="prefix")
parser.add_argument("-l", "--length_distribution_file_prefix", action="store",
                    dest="len_distr_file",
                    help="Output file with lengths distibutions",
                    default="length_distribution")

args = parser.parse_args()

examiner = GFFExaminer()

with open(args.gff, "r") as in_fd:
    pprint.pprint(examiner.parent_child_map(in_fd))

with open(args.gff, "r") as in_fd:
    record_dict = dict([(record.id, record) for record in GFF.parse(in_fd)])

gene_dict = OrderedDict({})
for record_id in record_dict:
    for feature in record_dict[record_id].features:
        if feature.type == "gene":
            gene_dict[feature.qualifiers["Name"][0]] = OrderedDict({})
            for sub_feature in feature.sub_features:
                gene_dict[feature.qualifiers["Name"][0]][sub_feature.type] = len(sub_feature)
        if feature.type in ("snoRNA", "ncRNA", "snRNA"):
            gene_dict[feature.qualifiers["Name"][0]] = OrderedDict({"ncRNA": len(feature)})

with open("%s_test.t" % args.prefix, "w") as out_fd:
    for gene in gene_dict:
        for sub_feature in gene_dict[gene]:
            out_fd.write("%s\t%s\t%i\n" % (gene, sub_feature, gene_dict[gene][sub_feature]))

lengths_dict = get_feature_lengths(record_dict)
count_dict = TwoLvlDict({})
for record in lengths_dict:
    count_dict[record] = {}
    for feature_type in lengths_dict[record]:
        count_dict[record][feature_type] = len(lengths_dict[record][feature_type])

count_dict.write("%s_counts.t" % args.prefix)
total_lengths = get_total_feature_lengths(lengths_dict, out_filename="%s_feature_lengths.t" % args.prefix)

white_list = ["five_prime_UTR", "three_prime_UTR", "CDS", "ncRNA"]
collapsed_dict = feature_lengths_collapse_records(lengths_dict,
                                                  synonym_dict={"snoRNA": "ncRNA", "snRNA": "ncRNA"})

for feature in collapsed_dict:
    collapsed_dict[feature] = np.array(collapsed_dict[feature])

bin_dict = {"five_prime_UTR": np.linspace(0, 900, 91), "three_prime_UTR": np.linspace(0, 1600, 81),
            "CDS": np.linspace(0, 16000, 81), "ncRNA": 40}


plt.figure(1, dpi=150, figsize=(24, 12))
index = 1
for feature_type in white_list:
    if feature_type not in collapsed_dict:
        continue
    plt.subplot(2, 2, index)
    plt.title(feature_type + " (%i)" % len(collapsed_dict[feature_type]))
    percentile_list = np.percentile(collapsed_dict[feature_type], [1, 5, 50, 95, 99])
    plt.hist(collapsed_dict[feature_type], bins=bin_dict[feature_type],
             label="Min: %i\nMax: %i\n1th percentile %i\n5th percentile %i\n50th percentile %i\n95th percentile %i\n99th percentile %i" %
                   (min(collapsed_dict[feature_type]),
                    max(collapsed_dict[feature_type]),
                    percentile_list[0], percentile_list[1], percentile_list[2], percentile_list[3], percentile_list[4]))
    plt.xlabel("Length")
    plt.ylabel("N")
    plt.legend()
    index += 1
plt.suptitle("Feature length distribution")
plt.savefig("%s_length_distribution.svg" % args.prefix)
plt.savefig("%s_length_distribution.eps" % args.prefix)
plt.close()