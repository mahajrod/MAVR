#!/usr/bin/env python
__author__ = 'mahajrod'

import os
from collections import OrderedDict
from Parsers.VCF import CollectionVCF
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

try:
    from BCBio import GFF
except:
    print("Please install bcbio-gff package to run this script. Exiting...")
    exit(0)


def variants_start_end(collection, left, right, record_dict, min_five_utr_len=10, skip_nonintergenic_variants=False):

    pre_UTR_positions = []
    UTR_positions = []
    CDS_positions = []

    variants_list = []
    for record_id in record_dict:
        for feature in record_dict[record_id].features:
            feature_strand = feature.location.strand
            if feature.type != "gene":
                continue

            for sub_feature in feature.sub_features:
                if sub_feature.type == "five_prime_UTR" and len(sub_feature) >= min_five_utr_len:
                    break
            else:
                continue
            #print(feature.sub_features)
            for sub_feature in feature.sub_features:
                strand = sub_feature.strand
                if sub_feature.type == "five_prime_UTR":
                    five_UTR_start = sub_feature.location.start + 1 if strand == +1 else sub_feature.location.end
                    pre_UTR_start = five_UTR_start - left if strand == +1 else five_UTR_start + 1
                    pre_UTR_end = five_UTR_start - 1 if strand == +1 else five_UTR_start + left
                    #print(sub_feature)
                    for variant in collection:
                        if record_id != variant.chrom:
                            continue
                        """
                        if (feature_strand == 1 and variant.ref != "C") or (feature_strand == -1 and variant.ref != "G"): #select mutations in un_transcribed thread
                            continue
                        """
                        if (feature_strand == 1 and variant.ref != "G") or (feature_strand == -1 and variant.ref != "C"): #select mutations in transcribed thread
                            continue

                        if variant.pos - 1 in sub_feature:
                            #print(five_UTR_start)
                            print(str(sub_feature.location) + "\t" + str(variant) + "\tUTR")
                            relative_position = float((variant.pos - five_UTR_start) * strand) * 100 / len(sub_feature)
                            if relative_position < 0:
                                print(variant)
                                print(sub_feature)
                            UTR_positions.append(relative_position)
                        elif pre_UTR_start <= variant.pos <= pre_UTR_end:
                            if skip_nonintergenic_variants and variant.info_dict["Ftype"] != ["igc"]:
                                continue
                            print(str(sub_feature.location) + "\t" + str(variant) + "\tpre_UTR")
                            relative_position = (variant.pos - five_UTR_start) * strand
                            if relative_position > 0:
                                print(pre_UTR_start, pre_UTR_end, five_UTR_start)
                                print(variant)
                                print(sub_feature.location)
                            pre_UTR_positions.append(relative_position)

                    continue
                if sub_feature.type != "CDS":
                    continue
                strand = sub_feature.strand
                CDS_start = sub_feature.location.start + 1 if strand == +1 else sub_feature.location.end

                region_start_start = CDS_start if strand == +1 else CDS_start - right
                region_start_end = CDS_start + right if strand == +1 else CDS_start

                for variant in collection:
                    if record_id != variant.chrom:
                        continue
                    """
                    if (feature_strand == 1 and variant.ref != "C") or (feature_strand == -1 and variant.ref != "G"): #select mutations in un_transcribed thread
                            continue
                    """

                    if (feature_strand == 1 and variant.ref != "G") or (feature_strand == -1 and variant.ref != "C"): #select mutations in transcribed thread
                        continue

                    if region_start_start <= variant.pos <= region_start_end:
                        CDS_positions.append((variant.pos - CDS_start) * strand)
                        print(str(sub_feature.location) + "\t" + str(variant) + "\tCDS")
    return pre_UTR_positions, UTR_positions, CDS_positions

if __name__ == "__main__":
    workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/clusters/all/all/"

    sample_set_names_list = ["PmCDA1_3d",
                             #"HAP",
                             "PmCDA1_sub1_3d",
                             "PmCDA1_6d",
                             #"HAP_sub1",
                             "PmCDA1_sub1_6d",
                             #"A1_3d",
                             #"A1_6d",
                             #"A3G_3d",
                             #"AID_3d",
                             #"AID_6d"
                             ]

    annotations = "/home/mahajrod/Genetics/Projects/desaminases/data/LAN210_v0.10m/annotations/merged_annotations_Nagalakshmi_tranf_to_LAN210_v0.10m.gff3"
    with open(annotations, "r") as in_fd:
        record_dict = dict([(record.id, record) for record in GFF.parse(in_fd)])
    pre_UTR_positions = OrderedDict({})
    UTR_positions = OrderedDict({})
    length_dict = OrderedDict({})
    CDS_positions = OrderedDict({})

    pre_UTR_hist_dict = OrderedDict({})
    UTR_hist_dict = OrderedDict({})
    CDS_hist_dict = OrderedDict({})
    os.chdir(workdir)

    left = 300
    right = 300
    bin_width = 5
    pre_UTR_bins = left / bin_width
    CDS_bins = np.linspace(0, right, right / bin_width + 1)
    UTR_bins = 10
    normed = True
    max_start = 0
    max_end = 0
    skip_nonintergenic_variants = True
    for sample_set in sample_set_names_list:
        print("Handling %s" % sample_set)
        vcf_file = "%s_good.vcf" % sample_set
        #start_hist_prefix = "%s_start_hist_r_%i_l_%i" % (sample_set, right, left)
        #end_hist_prefix = "%s_end_hist_r_%i_l_%i" % (sample_set, right, left)
        #gene_variants = "%s_gene_variants_r_%i_l_%i.t" % (sample_set, right, left)
        #variants, minus_variants = CollectionVCF(from_file=True, vcf_file=vcf_file).filter_by_expression("record.ref == 'C'") # C -> T variants
        #variants, minus_variants = CollectionVCF(from_file=True, vcf_file=vcf_file).filter_by_expression("(record.ref == 'C' and record.info_dict['Fstrand'][0] == 'P') or (record.ref == 'G' and record.info_dict['Fstrand'][0] == 'M')") # nontranscribed thread
        variants= CollectionVCF(from_file=True, vcf_file=vcf_file)
        pre_UTR_positions[sample_set], UTR_positions[sample_set], CDS_positions[sample_set] = \
            variants_start_end(variants, left, right, record_dict, min_five_utr_len=10,
                               skip_nonintergenic_variants=skip_nonintergenic_variants)
        length_dict[sample_set] = len(variants)
        #print(start_dict[sample_set])
        pre_UTR_hist_dict[sample_set] = list(np.histogram(pre_UTR_positions[sample_set], bins=pre_UTR_bins))
        UTR_hist_dict[sample_set] = list(np.histogram(UTR_positions[sample_set], bins=UTR_bins))
        CDS_hist_dict[sample_set] = list(np.histogram(CDS_positions[sample_set], bins=CDS_bins))
        print("UTR")
        print(UTR_positions[sample_set])
        print("blablabla")
        print(pre_UTR_hist_dict[sample_set][0])
        print(UTR_hist_dict[sample_set][0])
        print(CDS_hist_dict[sample_set][0])
        if normed:
            pre_UTR_hist_dict[sample_set][0] = pre_UTR_hist_dict[sample_set][0].astype(np.float32, copy=False)
            UTR_hist_dict[sample_set][0] = UTR_hist_dict[sample_set][0].astype(np.float32, copy=False)
            CDS_hist_dict[sample_set][0] = CDS_hist_dict[sample_set][0].astype(np.float32, copy=False)

            pre_UTR_hist_dict[sample_set][0] = pre_UTR_hist_dict[sample_set][0] / length_dict[sample_set]
            UTR_hist_dict[sample_set][0] = UTR_hist_dict[sample_set][0] / length_dict[sample_set]
            CDS_hist_dict[sample_set][0] = CDS_hist_dict[sample_set][0] / length_dict[sample_set]

        print("Normed")
        print(pre_UTR_hist_dict[sample_set][0])
        print(UTR_hist_dict[sample_set][0])
        print(CDS_hist_dict[sample_set][0])
        max_start = max(max_start, np.amax(pre_UTR_hist_dict[sample_set][0]),
                        np.amax(UTR_hist_dict[sample_set][0]),
                        np.amax(CDS_hist_dict[sample_set][0]) )
    plt.figure(1, dpi=300, figsize=(24, 8*len(sample_set_names_list)))

    max_start = 0.0225

    index = 0
    for sample_set in sample_set_names_list:
        plt.subplot(len(sample_set_names_list), 3, index * 3 + 1)
        plt.bar(pre_UTR_hist_dict[sample_set][1][:-1], pre_UTR_hist_dict[sample_set][0], width=bin_width)
        plt.xlim(xmin=-left, xmax=0)
        plt.ylim(ymax=max_start)
        plt.axhline(0.02, color='y')
        plt.axhline(0.01, color='k')
        plt.axhline(0.005, color='r')
        plt.axhline(0.0025, color='g')
        plt.title(sample_set + " pre 5' UTR")

        plt.subplot(len(sample_set_names_list), 3, index * 3 + 2)
        plt.bar(UTR_hist_dict[sample_set][1][:-1], UTR_hist_dict[sample_set][0], width=10, color="#6686ED")
        plt.xlim(xmin=0, xmax=100)
        plt.ylim(ymax=max_start)
        plt.axhline(0.02, color='y')
        plt.axhline(0.01, color='k')
        plt.axhline(0.005, color='r')
        plt.axhline(0.0025, color='g')
        plt.title(sample_set + " 5' UTR (% pos)")

        plt.subplot(len(sample_set_names_list), 3, index * 3 + 3)
        plt.bar(CDS_hist_dict[sample_set][1][:-1], CDS_hist_dict[sample_set][0], width=bin_width)
        plt.xlim(xmin=0, xmax=right)
        plt.ylim(ymax=max_start)
        plt.axhline(0.02, color='y')
        plt.axhline(0.01, color='k')
        plt.axhline(0.005, color='r')
        plt.axhline(0.0025, color='g')
        plt.title(sample_set + " CDS")
        index += 1

    suffix = "pre_five_UTR_only_intergenic" if skip_nonintergenic_variants else "all"
    plt.savefig("TSS_UTR_CDS_start_all_r_%i_l_%i_bin_width_%i_%s_transcribed_variants.svg" % (right, left, bin_width, suffix))
    plt.savefig("TSS_UTR_CDS_start_all_r_%i_l_%i_bin_width_%i_%s_transcribed_variants.eps" % (right, left, bin_width, suffix))
    plt.close()