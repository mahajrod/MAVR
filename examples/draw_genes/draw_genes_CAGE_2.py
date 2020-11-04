#!/usr/bin/env python

from collections import OrderedDict
from Bio import SeqIO
from BCBio import GFF
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as tck

from RouToolPa.GeneralRoutines.File import read_tsv_as_columns_dict


def check_location_intersection(location1, location2):
    if location1.start in location2 or location1.end in location2:
        return True
    if location2.start in location1 or location2.end in location1:
        return True
    return False

CAGE_dir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/CAGE/GSE66284/"
#annotation_file = "/home/mahajrod/Reference_genomes/drosophila_melanogaster/r6.03/gtf/sbr.gtf"
annotation_file = "/home/mahajrod/Reference_genomes/drosophila_melanogaster/r6.03/gtf/sbr_locus.gtf"
description_file = CAGE_dir + "good_samples_description.t"
length_file = "length_sam_file.len"
histo_file = "nxf_region_filtered_sorted_filtered.histo"
with open(annotation_file, "r") as ann_fd:
    annotations_dict = SeqIO.to_dict(GFF.parse(ann_fd))


samples_description_dict = read_tsv_as_columns_dict(description_file)
number_of_samples = len(samples_description_dict["sample"])
samples_sizes_dict = OrderedDict()
sample_list = ["SRR488284",
               "SRR488285",
               "SRR488308",
               "SRR488309"]
sample_size = {
               "SRR488284": 7431086,
               "SRR488285": 7900633,
               "SRR488308": 7638501,
               "SRR488309": 8406294}
samples_dict = OrderedDict({})
rev_samples_dict = OrderedDict({})
print(samples_description_dict.keys())
for sample in samples_description_dict["sample"]:
    with open(CAGE_dir + sample + "/" + length_file, "r") as len_fd:
        samples_sizes_dict[sample] = int(len_fd.readline().strip().split()[0]) - 10
    samples_dict[sample] = [[], []]
    rev_samples_dict[sample] = [[], []]
    with open(CAGE_dir + sample + "/" + histo_file, "r") as in_fd:
        in_fd.readline()
        for line in in_fd:
            line_list = line.strip().split()
            if line_list[0] == "-1":
                samples_dict[sample][0].append(int(line_list[1]))
                samples_dict[sample][1].append(float(line_list[2]) * float(1000000) / float(samples_sizes_dict[sample]))
            else:
                rev_samples_dict[sample][0].append(int(line_list[1]))
                rev_samples_dict[sample][1].append(float(line_list[2]) * float(1000000) / float(samples_sizes_dict[sample]))



five_prime_utr_type = "5UTR"
three_prime_utr_type = "3UTR"
exon_type = "exon"
cds_type = "CDS"


cds_height = 14
utr_height = 8
print (samples_dict)
vertical_distance = 14

borders = (10832700, 10847386)

ymax = max([max(samples_dict[sample][1]) for sample in samples_dict])

figure = plt.figure(1, (15, 2.5 * (number_of_samples + 1)), dpi=200)
subplot_index = 1
subplot_dict = OrderedDict()
ymax_dict = OrderedDict()
axes_dict = OrderedDict()

for chrom in annotations_dict:
        for feature in annotations_dict[chrom].features:
            if feature.qualifiers["ID"][0] != "FBtr0073411":
                continue
            five_prime_utr_location_list = []
            three_prime_utr_location_list = []
            exon_location_list = []
            cds_location_list = []
            print feature.qualifiers["ID"][0]
            #print feature.sub_features
            for sub_feature in feature.sub_features:
                if sub_feature.type == five_prime_utr_type:
                    five_prime_utr_location_list.append([sub_feature.location.start + 1, int(sub_feature.location.end)])
                elif sub_feature.type == three_prime_utr_type:
                    three_prime_utr_location_list.append([sub_feature.location.start + 1, int(sub_feature.location.end)])
                elif sub_feature.type == exon_type:
                    strand = sub_feature.location.strand
                    exon_location_list.append([sub_feature.location.start + 1, int(sub_feature.location.end)])
                elif sub_feature.type == cds_type:
                    cds_location_list.append([sub_feature.location.start + 1, int(sub_feature.location.end)])

            five_prime_utr_location_list.sort()
            exon_location_list.sort()
            three_prime_utr_location_list.sort()
            cds_location_list.sort()

print(three_prime_utr_location_list)

for sample in samples_dict:
    #print(number_of_samples)
    subplot = plt.subplot(number_of_samples + 1, 1, subplot_index)
    subplot_dict[sample] = subplot
    axes = figure.add_subplot(number_of_samples + 1, 1, subplot_index)
    axes_dict[sample] = axes
    #axes.spines['bottom'].set_color('none')
    axes.spines['right'].set_color('none')
    #axes.spines['left'].set_color('none')
    axes.spines['top'].set_color('none')
    plt.ylabel("TPM", fontsize=8)
    #axes.yaxis.set_ticks([0, 20, 40])
    plt.yticks(fontsize=7)
    plt.text(-0.1, 0.5, sample, rotation=90, fontweight="bold", transform=subplot.transAxes, fontsize=9,
                         horizontalalignment='center',
                         verticalalignment='center')
    ymax = max(max(samples_dict[sample][1]), max(rev_samples_dict[sample][1]) if rev_samples_dict[sample][1] else 0)
    ymax_dict[sample] = ymax

    axes.xaxis.set_ticks([])

    axes.yaxis.tick_left()

    for cds_start, cds_end in cds_location_list:
        axes.add_artist(plt.Rectangle((cds_start, 0), cds_end - cds_start + 1, ymax_dict[sample] + 1,
                                              facecolor="black", linewidth=0, alpha=0.35))
    for three_utr_start, three_utr_end in three_prime_utr_location_list:
        axes.add_artist(plt.Rectangle((three_utr_start, 0), three_utr_end - three_utr_start + 1,
                                              ymax_dict[sample] + 1, facecolor="#cacaca", linewidth=0, alpha=0.35))

    for five_utr_start, five_utr_end in five_prime_utr_location_list:
        axes.add_artist(plt.Rectangle((five_utr_start, 0), five_utr_end - five_utr_start + 1,
                                              ymax_dict[sample] + 1, facecolor="#606060", linewidth=0, alpha=0.35))

    plt.bar(samples_dict[sample][0], samples_dict[sample][1], width=1, color="red", edgecolor="red")
    plt.bar(rev_samples_dict[sample][0], rev_samples_dict[sample][1], width=1, color="green", edgecolor="green")

    ax = plt.gca()

    plt.xlim(xmin=borders[0], xmax=borders[1])
    plt.ylim(ymin=0, ymax=ymax + 1)
    ax.invert_xaxis()
    subplot_index += 1

subplot = plt.subplot(number_of_samples + 1, 1, subplot_index)

x_formatter = tck.ScalarFormatter(useOffset=False)
x_formatter.set_scientific(False)

axes = figure.add_subplot(number_of_samples + 1, 1, subplot_index)

axes.get_yaxis().set_visible(False)
axes.get_xaxis().set_visible(False)
#axes.xaxis.set_major_formatter(x_formatter)

axes.spines['bottom'].set_color('none')
axes.spines['right'].set_color('none')
axes.spines['left'].set_color('none')
axes.spines['top'].set_color('none')
"""
axes.xaxis.set_ticks_position('bottom')
axes.xaxis.set_major_formatter(tck.FuncFormatter(lambda x, p: format(int(x), ',')))
"""
left_coord = None
right_coord = None

locations_dict = {1: [], -1: []}

for chrom in annotations_dict:
        for feature in annotations_dict[chrom].features:
            five_prime_utr_location_list = []
            three_prime_utr_location_list = []
            exon_location_list = []
            cds_location_list = []
            print feature.qualifiers["ID"][0]
            #print feature.sub_features
            for sub_feature in feature.sub_features:
                if sub_feature.type == five_prime_utr_type:
                    five_prime_utr_location_list.append([sub_feature.location.start + 1, int(sub_feature.location.end)])
                elif sub_feature.type == three_prime_utr_type:
                    three_prime_utr_location_list.append([sub_feature.location.start + 1, int(sub_feature.location.end)])
                elif sub_feature.type == exon_type:
                    strand = sub_feature.location.strand
                    exon_location_list.append([sub_feature.location.start + 1, int(sub_feature.location.end)])
                elif sub_feature.type == cds_type:
                    cds_location_list.append([sub_feature.location.start + 1, int(sub_feature.location.end)])

            five_prime_utr_location_list.sort()
            exon_location_list.sort()
            three_prime_utr_location_list.sort()
            cds_location_list.sort()

            feature_location = feature.location
            # find line in picture
            if locations_dict[strand]:
                for i in range(0, len(locations_dict[strand])):
                    for location in locations_dict[strand][i]:
                        if check_location_intersection(feature_location, location):
                            break
                    else:
                        line = i
                        locations_dict[strand][i].append(feature_location)
                        break
                else:
                    line = i + 1
                    locations_dict[strand].append([feature_location])
            else:
                locations_dict[strand].append([feature_location])
                line = 0

            left_coord = exon_location_list[0][0] if left_coord is None else min(left_coord, exon_location_list[0][0])
            right_coord = exon_location_list[-1][-1] if right_coord is None else max(right_coord, exon_location_list[-1][-1])

            vertical_shift = (vertical_distance + cds_height) * strand * line + 2 * cds_height

            cds_vertical_pos = - 2*vertical_distance - cds_height if strand == -1 else vertical_distance
            utr_vertical_pos = - 2*vertical_distance - (cds_height+utr_height)/2 if strand == -1 else vertical_distance + (cds_height-utr_height)/2
            cds_vertical_pos += vertical_shift
            utr_vertical_pos += vertical_shift

            intron_line_y = cds_vertical_pos + cds_height/2


            """
            for exon_start, exon_end in exon_location_list:
                axes.add_patch(plt.Rectangle((exon_start, vertical_pos), exon_end - exon_start + 1, height, facecolor=None, edgecolor='black'))
            """
            for i in range(0, len(exon_location_list)-1):
                axes.add_line(plt.Line2D([exon_location_list[i][-1] + 1, exon_location_list[i+1][0] - 1],
                                         [intron_line_y, intron_line_y], linewidth=1, color='black'))
            for cds_start, cds_end in cds_location_list:
                axes.add_artist(plt.Rectangle((cds_start, cds_vertical_pos), cds_end - cds_start + 1, cds_height,
                                              facecolor="black", linewidth=0))

            for three_utr_start, three_utr_end in three_prime_utr_location_list:
                axes.add_artist(plt.Rectangle((three_utr_start, utr_vertical_pos), three_utr_end - three_utr_start + 1,
                                              utr_height, facecolor="#cacaca", linewidth=0))

            for five_utr_start, five_utr_end in five_prime_utr_location_list:
                axes.add_artist(plt.Rectangle((five_utr_start, utr_vertical_pos), five_utr_end - five_utr_start + 1,
                                              utr_height, facecolor="#606060", linewidth=0))

            name_x = exon_location_list[0][0]
            #name_x = five_prime_utr_location_list[0][0] if strand == 1 else five_prime_utr_location_list[-1][0]
            name_y = cds_vertical_pos + cds_height * 1.25

            #axes.text(name_x, name_y, feature.sub_features[0].qualifiers["gene_symbol"][0], style='italic', fontsize=9)

axes.xaxis.set_ticks([])
axes.invert_xaxis()
dps = figure.dpi_scale_trans.inverted()
bbox = axes.get_window_extent().transformed(dps)
width, height = bbox.width, bbox.height
xmin = left_coord - 10
xmax = right_coord + 10

plus_lines = len(locations_dict[1])
minus_lines = len(locations_dict[-1])

ymin = - (minus_lines * (vertical_distance + cds_height) + 2 * cds_height)
ymax = plus_lines * (vertical_distance + cds_height) + 2 * cds_height

#axes.text(xmin, cds_height/2, "chr %s" % annotations_dict.keys()[0], fontsize=10, fontweight='bold')
#plt.xticks(fontsize=10)

ax = plt.gca()

plt.xlim(xmin=borders[0], xmax=borders[1])
plt.ylim(ymin=ymin, ymax=ymax)
ax.invert_xaxis()
for extension in ".jpg", ".svg", ".png":
    plt.savefig("CAGE_sbr_locus_all_TPM_both_strands%s" % extension)
