#!/usr/bin/env python
__author__ = 'mahajrod'

import os
import numpy as np

from Parsers.VCF import CollectionVCF
from Parsers.CCF import CollectionCCF
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams
from BCBio import GFF


def location_pie(self, annotation_colors=[],
                     ref_genome=None, explode=True, annotation_black_list=[],
                     allow_several_counts_of_record=False,
                     counts_filename="location_counts.t",
                     counts_dir="location_counts",
                     legend_font=7,
                     combine_mixed=False):
        reference_colors = {"CDS": "#FBFD2B",    # yellow
                            "5'_UTR": "#FF000F",
                            "five_prime_UTR": "#FF000F",     # red
                            "3'_UTR": "#000FFF",
                            "three_prime_UTR": "#000FFF",     # blue
                            "igc": "#4ED53F",     # green
                            "ncRNA": 'cyan',
                            "other": "#ADB696"
                            }

        if annotation_colors:
            reference_colors = annotation_colors
            reference_colors["other"] = "#ADB696"

        count_locations_dict = self.count_locations(annotation_black_list=annotation_black_list,
                                                    allow_several_counts_of_record=allow_several_counts_of_record,
                                                    out_filename=counts_filename,
                                                    write=True,
                                                    count_dir=counts_dir)

        index = 1
        all_labels = []
        all_counts = []
        all_colors = []

        for region in count_locations_dict:
            labels = []
            counts = []
            colors = []
            #print(count_locations_dict)
            if combine_mixed:
                labels.append("other")
                counts.append(0)
                colors.append(reference_colors["other"])
            for label in count_locations_dict[region]:
                #print(count_locations_dict[region])

                if count_locations_dict[region][label] == 0 or label in annotation_black_list:
                    continue
                if combine_mixed and "/" in label:
                    counts[0] += count_locations_dict[region][label]
                    continue

                labels.append(label)
                counts.append(count_locations_dict[region][label])
                if label not in reference_colors:
                    colors.append(reference_colors["other"])
                else:
                    colors.append(reference_colors[label])

            for i in range(0, len(labels)):
                if labels[i] not in all_labels:
                    all_labels.append(labels[i])
                    all_counts.append(counts[i])
                    all_colors.append(colors[i])
                else:
                    label_index = all_labels.index(labels[i])
                    all_counts[label_index] += counts[i]


        all_explodes = np.zeros(len(all_counts))
        if explode and all_counts:
            max_count_index = all_counts.index(max(all_counts))
            all_explodes[max_count_index] = 0.1
        if len(all_labels) > 0:
            max_label_length = max([len(x) for x in all_labels])
            max_count = max(all_counts)
            max_letters = 1
            while int(max_count / 10**max_letters) != 0:
                max_letters += 1

        patches, texts = plt.pie(all_counts, explode=all_explodes, colors=all_colors,
                                            shadow=True, startangle=90, radius=2) #labels=all_labels,

        #patches, texts, autotexts = plt.pie(all_counts, explode=all_explodes, colors=all_colors,
        #                                    shadow=True, startangle=90, autopct='%1.1f%%', radius=4) #labels=all_labels

        #colors = ['yellowgreen','red','gold','lightskyblue','white','lightcoral','blue','pink', 'darkgreen','yellow','grey','violet','magenta','cyan']
        porcent = 100 * np.array(all_counts).astype(np.float32, copy=False)/sum(all_counts)

        labels = ['{0}  -  {1} ({2:1.1f}%)'.format(i.ljust(max_label_length), str(j).ljust(max_letters), k) for i, j, k in zip(all_labels, all_counts, porcent)]

        if len(all_labels) > 0:
            sort_legend = True
            if sort_legend:
                patches, labels, dummy = zip(*sorted(zip(patches, labels, all_counts),
                                                     key=lambda x: x[2],
                                                     reverse=True))

        plt.legend(patches, labels, loc='center left',  fontsize=legend_font, bbox_to_anchor=(-0.25, 0.5)) # bbox_to_anchor=(-0.1, 1.),
            # Set aspect ratio to be equal so that pie is drawn as a circle.
        plt.axis('equal')
        #plt.savefig("%s/%s" % (plot_dir, full_genome_pie_filename), bbox_inches='tight')

if __name__ == "__main__":
    workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/clusters/all/all/"
    letter_list_part1 = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]
    os.chdir(workdir)
    sample_set_names_list = ["PmCDA1_3d",
                             "PmCDA1_sub1_3d",
                             "PmCDA1_6d",
                             "PmCDA1_sub1_6d",
                             "HAP",
                             "HAP_sub1"
                             ]
    annotation_black_list = ["gene", "region", "ARS", "long_terminal_repeat",
                             "noncoding_exon", "intron", "repeat_region", "telomere", "gene_cassette",
                              "five_prime_UTR_intron", "pseudogene", "transposable_element_gene", "LTR_retrotransposon"]
    gff_file = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/annotations/merged_annotations_Nagalakshmi_tranf_to_LAN210_v0.10m.gff3"
    annotations_dict = {}
    annotation_synonym_dict = {"three_prime_UTR": "3'_UTR",
                               "five_prime_UTR": "5'_UTR",
                               "snoRNA": "ncRNA",
                               "snRNA": "ncRNA"
                               }
    with open(gff_file) as gff_fd:
        for record in GFF.parse(gff_fd):
            annotations_dict[record.id] = record


    rcParams.update({'font.size': 7})
    plt.figure(1, dpi=300, figsize=(8, 6))
    index = 1
    for sample, letter in zip(sample_set_names_list, letter_list_part1):
        collection = CollectionVCF(from_file=True, vcf_file=sample + "_good.vcf")

        plt.subplot(3, 2, index)
        collection.get_location(annotations_dict, use_synonym=True, synonym_dict=annotation_synonym_dict)
        location_pie(collection, annotation_colors=[],
                     ref_genome=None, explode=True, annotation_black_list=annotation_black_list,
                     allow_several_counts_of_record=False,
                     counts_filename="location_counts.t",
                     counts_dir="location_counts",
                     legend_font=6,
                     combine_mixed=True
                     )
        plt.title("%s. %s" % (letter, sample), fontweight='bold')

        index += 1
    for format_ext in ["svg", "eps", "pdf", "png"]:
        plt.savefig("good_mutation_pie_mixed_combined.%s" % format_ext, bbox_inches='tight')
    plt.close()

