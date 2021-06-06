#!/usr/bin/env python

from collections import OrderedDict
from Bio import SeqIO
try:
    from BCBio import GFF
except:
    print("Please install bcbio-gff package to run this script. Exiting...")
    exit(0)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as tck

import mpl_toolkits.axisartist as axisartist
from matplotlib.patches import Rectangle


def check_location_intersection(location1, location2):
    if location1.start in location2 or location1.end in location2:
        return True
    if location2.start in location1 or location2.end in location1:
        return True
    return False

annotation_file = "/home/mahajrod/Reference_genomes/drosophila_melanogaster/r6.03/gtf/sbr.gtf"
with open(annotation_file, "r") as ann_fd:
    annotations_dict = SeqIO.to_dict(GFF.parse(ann_fd))

CAGE_dir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/CAGE/GSE66284/"
annotation_file = "/home/mahajrod/Reference_genomes/drosophila_melanogaster/r6.03/gtf/sbr_locus.gtf"
description_file = CAGE_dir + "good_samples_description.t"
length_file = "length_sam_file.len"
histo_file = "nxf_region_filtered_sorted_filtered.histo"

with open(annotation_file, "r") as ann_fd:
    annotations_dict = SeqIO.to_dict(GFF.parse(ann_fd))

samples_description_dict = {"sample": ["SRR488285",
                                       #"SRR488282",
                                       #"SRR488271",
                                       "SRR488272"],
                            "description": ["Adult Mated Male 4 days Post-eclosion Testes",
                                            #"Adult Mated Female 4 days Post-eclosion Ovaries",
                                            #"Adult Mated Female 1 day Post-eclosion Heads",
                                            "Adult Mated Male 1 day Post-eclosion Heads"]
                            }

description_dict = {"SRR488285": "Testes\n(4 days male)",
                    #"SRR488282": "Ovaries\n(4 days female)",
                    #"SRR488271": "Heads\n(1 day female)",
                    "SRR488272": "Heads\n(1 day male)"
                    }

number_of_samples = len(samples_description_dict["sample"])
samples_sizes_dict = OrderedDict()

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


cds_height = 28
utr_height = 16
print (samples_dict)
vertical_distance = 14

borders = [(10846450, 10847350), (10836650, 10837550)]

ymax = max([max(samples_dict[sample][1]) for sample in samples_dict])

figure = plt.figure(1, (5, 5), dpi=300)
subplot_index = 1
gene_start = 10847092
cage_peaks_dict = {"SRR488272": [10846888, 10847111],
                   "SRR488285": [10837480, 10846888, 10847093]}


def calc_relative_coordinate_label(peak_coord):
    return str(gene_start - peak_coord + 1) if gene_start - peak_coord >= 0 else str(gene_start - peak_coord)

for sample in samples_dict:
    for region_index in [0, 1]:
        subplot = plt.subplot(number_of_samples + 1, 2, subplot_index)
        axes = figure.add_subplot(number_of_samples + 1, 2, subplot_index)
        #axes.spines['bottom'].set_color('none')
        axes.spines['right'].set_color('none')
        #axes.spines['left'].set_color('none')
        axes.spines['top'].set_color('none')
        if region_index == 1:
            axes.spines['left'].set_color('none')
            axes.yaxis.set_ticks([])
        if region_index == 0:
            plt.ylabel("TPM", fontsize=9)
            #axes.yaxis.set_ticks([0, 20, 40])
            plt.yticks(fontsize=9)
            plt.text(-0.30, 0.5, description_dict[sample], rotation=90, fontweight="bold",
                     transform=subplot.transAxes, fontsize=10,
                     horizontalalignment='center',
                     verticalalignment='center')
        plt.bar(samples_dict[sample][0], samples_dict[sample][1], width=1)
        plt.xlim(xmin=borders[region_index][0], xmax=borders[region_index][1])
        plt.ylim(ymin=0, ymax=ymax)
        axes.xaxis.set_ticks_position('bottom')
        axes.xaxis.set_ticks(cage_peaks_dict[sample])
        axes.xaxis.set_ticklabels(map(calc_relative_coordinate_label, cage_peaks_dict[sample]), fontsize=9)
        axes.xaxis.set_tick_params(direction='out')
        axes.yaxis.tick_left()
        plt.xlim(xmin=borders[region_index][0], xmax=borders[region_index][1])
        plt.ylim(ymin=0, ymax=ymax)
        ax = plt.gca()
        ax.invert_xaxis()
        subplot_index += 1
for subplot_index, limits in zip([9, 10], borders):
    subplot = plt.subplot((number_of_samples + 1)*2, 2, subplot_index)

    x_formatter = tck.ScalarFormatter(useOffset=False)
    x_formatter.set_scientific(False)

    axes = figure.add_subplot((number_of_samples + 1)*2, 2, subplot_index)

    axes.get_yaxis().set_visible(False)
    axes.get_xaxis().set_visible(True)
    axes.xaxis.set_ticks_position('bottom')
    #axes.xaxis.set_major_formatter(x_formatter)

    #axes.spines['bottom'].set_color('none')
    axes.spines['right'].set_color('none')
    axes.spines['left'].set_color('none')
    axes.spines['top'].set_color('none')
    """

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
            print(feature.qualifiers["ID"][0])
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

            vertical_shift = 50 + (vertical_distance + cds_height) * strand * line + 2 * cds_height

            cds_vertical_pos = - 2*vertical_distance - cds_height if strand == -1 else vertical_distance
            utr_vertical_pos = - 2*vertical_distance - (cds_height+utr_height)/2 if strand == -1 else vertical_distance + (cds_height-utr_height)/2
            #cds_vertical_pos += vertical_shift
            #utr_vertical_pos += vertical_shift

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



    #axes.xaxis.set_ticks([])
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
    """
    # manual arrowhead width and length
    hw = 1./20.*(ymax-ymin)
    hl = 1./20.*(xmax-xmin)
    lw = 1.  # axis line width
    ohg = 0.3  # arrow overhang

    # compute matching arrowhead length and width
    yhw = hw / (ymax-ymin) * (xmax-xmin) * height / width
    yhl = hl / (xmax-xmin) * (ymax-ymin) * width / height

    # draw x and y axis

    axes.arrow(xmin, 0, xmax-xmin, 0., fc='k', ec='k', lw=lw,
               head_width=hw, head_length=hl, overhang=ohg,
               length_includes_head=True, clip_on=False)
    """
    #axes.text(xmin, cds_height/2, "chr %s" % annotations_dict.keys()[0], fontsize=10, fontweight='bold')
    #plt.xticks(fontsize=10)


    ax = plt.gca()
    ax.set_xticks([10836693, 10837214, 10846492, 10847092])
    ax.set_xticklabels(["10400", "9879", "601", "1"], fontsize=9)
    ax.xaxis.set_tick_params(direction='out')
    #ax.xaxis.set_ticks([10847092], [1])


    if subplot_index == 9:
        plt.text(0.6, 0.15, "Exon 1", rotation=0, fontweight="bold", transform=subplot.transAxes, fontsize=9,
                         horizontalalignment='center',
                         verticalalignment='center')
    if subplot_index == 10:
        plt.text(0.65, 0.15, "Exon 4", rotation=0, fontweight="bold", transform=subplot.transAxes, fontsize=9,
                         horizontalalignment='center',
                         verticalalignment='center')
    plt.xlim(xmin=limits[0], xmax=limits[1])
    plt.ylim(ymin=ymin, ymax=ymax)
    ax.invert_xaxis()
plt.subplots_adjust(hspace=0.25, wspace=0.10, top=0.95, left=0.17, right=0.97, bottom=0.01)
for extension in ".jpg", ".svg", ".png":
    plt.savefig("CAGE_sbr_locus_TPM_article_2_samples_same_scale%s" % extension)
