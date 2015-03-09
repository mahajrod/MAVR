#!/usr/bin/env python

from Bio import SeqIO
from BCBio import GFF

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

annotation_file = "/home/mahajrod/Reference_genomes/drosophila_melanogaster/r6.03/gtf/sbr_locus.gtf"
with open(annotation_file, "r") as ann_fd:
    annotations_dict = SeqIO.to_dict(GFF.parse(ann_fd))


five_prime_utr_type = "5UTR"
three_prime_utr_type = "3UTR"
exon_type = "exon"
cds_type = "CDS"


cds_height = 14
utr_height = 8

vertical_distance = 14

figure = plt.figure(1, (10, 3), dpi=300)
subplot = plt.subplot(1, 1, 1)

x_formatter = tck.ScalarFormatter(useOffset=False)
x_formatter.set_scientific(False)

axes = figure.add_subplot(111)
axes.get_yaxis().set_visible(False)
axes.xaxis.set_major_formatter(x_formatter)
axes.spines['bottom'].set_position('zero')
axes.spines['right'].set_color('none')
axes.spines['left'].set_color('none')
axes.spines['top'].set_color('none')
axes.xaxis.set_ticks_position('bottom')
axes.xaxis.set_major_formatter(tck.FuncFormatter(lambda x, p: format(int(x), ',')))

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

        vertical_shift = (vertical_distance + cds_height) * strand * line

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

        axes.text(name_x, name_y, feature.sub_features[0].qualifiers["gene_symbol"][0], style='italic', fontsize=9)





plt.title(r"Structure of $sbr$ locus")
axes.xaxis.set_ticks([left_coord, right_coord])
dps = figure.dpi_scale_trans.inverted()
bbox = axes.get_window_extent().transformed(dps)
width, height = bbox.width, bbox.height
xmin = left_coord - 1000
xmax = right_coord + 1000

plus_lines = len(locations_dict[1])
minus_lines = len(locations_dict[-1])

ymin = - (minus_lines * (vertical_distance + cds_height) + 2 * cds_height)
ymax = plus_lines * (vertical_distance + cds_height) + 2 * cds_height
print (plus_lines, minus_lines)
print(ymin, ymax)

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

axes.text(xmin, cds_height/2, "chr %s" % annotations_dict.keys()[0], fontsize=10, fontweight='bold')
plt.xticks(fontsize=10)
plt.xlim(xmin=xmin, xmax=xmax)
plt.ylim(ymin=ymin, ymax=ymax)
for extension in ".jpg", ".svg", ".png":
    plt.savefig("sbr_locus%s" % extension)
