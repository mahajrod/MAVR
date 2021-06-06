#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox, TransformedBbox, \
     blended_transform_factory

from mpl_toolkits.axes_grid1.inset_locator import BboxPatch, BboxConnector,\
     BboxConnectorPatch

import os
from collections import OrderedDict
import numpy as np

from Parsers.VCF import ReferenceGenome, CollectionVCF, ref_alt_variants

from Parsers.GFF import CollectionGFF
try:
    from BCBio import GFF
except:
    print("Please install bcbio-gff package to run this script. Exiting...")
    exit(0)

def connect_bbox(bbox1, bbox2,
                 loc1a, loc2a, loc1b, loc2b,
                 prop_lines, prop_patches=None):
    if prop_patches is None:
        prop_patches = prop_lines.copy()
        prop_patches["alpha"] = prop_patches.get("alpha", 1)*0.2

    c1 = BboxConnector(bbox1, bbox2, loc1=loc1a, loc2=loc2a, **prop_lines)
    c1.set_clip_on(False)
    c2 = BboxConnector(bbox1, bbox2, loc1=loc1b, loc2=loc2b, **prop_lines)
    c2.set_clip_on(False)

    bbox_patch1 = BboxPatch(bbox1, **prop_patches)
    bbox_patch2 = BboxPatch(bbox2, **prop_patches)

    p = BboxConnectorPatch(bbox1, bbox2,
                           #loc1a=3, loc2a=2, loc1b=4, loc2b=1,
                           loc1a=loc1a, loc2a=loc2a, loc1b=loc1b, loc2b=loc2b,
                           **prop_patches)
    p.set_clip_on(False)

    return c1, c2, bbox_patch1, bbox_patch2, p


def zoom_effect01(ax1, ax2, xmin, xmax, **kwargs):
    """
    ax1 : the main axes
    ax1 : the zoomed axes
    (xmin,xmax) : the limits of the colored area in both plot axes.

    connect ax1 & ax2. The x-range of (xmin, xmax) in both axes will
    be marked.  The keywords parameters will be used ti create
    patches.

    """

    trans1 = blended_transform_factory(ax1.transData, ax1.transAxes)
    trans2 = blended_transform_factory(ax2.transData, ax2.transAxes)

    bbox = Bbox.from_extents(xmin, 0, xmax, 1)

    mybbox1 = TransformedBbox(bbox, trans1)
    mybbox2 = TransformedBbox(bbox, trans2)

    prop_patches=kwargs.copy()
    prop_patches["ec"]="none"
    prop_patches["alpha"]=0.2

    c1, c2, bbox_patch1, bbox_patch2, p = \
        connect_bbox(mybbox1, mybbox2,
                     loc1a=3, loc2a=2, loc1b=4, loc2b=1,
                     prop_lines=kwargs, prop_patches=prop_patches)

    ax1.add_patch(bbox_patch1)
    ax2.add_patch(bbox_patch2)
    ax2.add_patch(c1)
    ax2.add_patch(c2)
    ax2.add_patch(p)

    return c1, c2, bbox_patch1, bbox_patch2, p


def rainfall_plot(self, region, base_colors=[], facecolor="#D6D6D6",
                  ref_genome=None, masked_regions=None, min_gap_length=10, draw_gaps=False):
        print("Drawing rainfall plot...")
        plot_dir = "rainfall_plot"
        reference_colors = {"A": "#FBFD2B",    # yellow
                            "C": "#FF000F",     # red
                            "G": "#000FFF",     # blue
                            "T": "#4ED53F",     # green
                            "INDEL": "#000000"  # black
                            }
        if base_colors:
            reference_colors = base_colors


        regions_dict = self._split_regions()
        num_of_regions = len(regions_dict)
        positions_dict = {}
        distances_dict = {}
        region_reference_dict = {}


        positions_dict[region] = np.array([record.pos for record in regions_dict[region]])
        #np.ediff1d return differences between consecutive elements in array, then 0 is added to the beginning
        distances_dict[region] = np.insert(np.ediff1d(positions_dict[region]), 0, 0)
        maximum = np.amax(distances_dict[region])
        region_reference_dict[region] = OrderedDict({"A": [[], []],
                                                         "C": [[], []],
                                                         "G": [[], []],
                                                         "T": [[], []],
                                                         "INDEL": [[], []]})
        for i in range(0, len(regions_dict[region])):
            region_reference_dict[region][self._reference(regions_dict[region][i])][0].append(positions_dict[region][i])
            region_reference_dict[region][self._reference(regions_dict[region][i])][1].append(distances_dict[region][i])
            if draw_gaps:
                if ref_genome:
                    for gap in ref_genome.gaps_dict[region]:
                        plt.gca().add_patch(plt.Rectangle((gap.location.start, 1),
                                                              gap.location.end - gap.location.start,
                                                              1024*32, facecolor="#777777", edgecolor='none'))
                # masked regions should be SeqRecord dict
                if masked_regions:
                    for feature in masked_regions[region].features:
                         plt.gca().add_patch(plt.Rectangle((int(feature.location.start)+1, 1),
                                                              feature.location.end - feature.location.start,
                                                              1024*32, facecolor="#aaaaaa", edgecolor='none'))

            for reference in region_reference_dict[region]:
                plt.plot(region_reference_dict[region][reference][0],
                             region_reference_dict[region][reference][1],
                             color=reference_colors[reference],
                             marker='.', linestyle='None', label=reference)

            #plt.title("Region %s" % region)
            plt.xlabel("Coordinate")
            plt.ylabel("Distanse")
            #plt.ylim(ymin=0)
            plt.axhline(y=100, color="#000000")
            plt.axhline(y=1000, color="#000000")
            plt.axhline(y=500, color="purple")
            plt.axhline(y=10, color="#000000")

                #if ref_genome:
                #    plt.xlim(xmax=len(ref_genome.reference_genome[region]))

        return region_reference_dict, maximum

if __name__ == "__main__":
    workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/clusters/all/all/"

    sample_set_names_list = ["PmCDA1_3d",
                             #"HAP",
                             #"PmCDA1_sub1_3d",
                             #"PmCDA1_6d",
                             #"HAP_sub1",
                             #"PmCDA1_sub1_6d"
                             ]
    reference = ReferenceGenome("/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/LAN210_v0.10m.fasta",
                                index_file="/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/LAN210_v0.10m.idx")

    reference_colors = {"A": "#FBFD2B",    # yellow
                            "C": "#FF000F",     # red
                            "G": "#000FFF",     # blue
                            "T": "#4ED53F",     # green
                            "INDEL": "#000000"  # black
                            }
    reference.find_gaps()
    os.chdir(workdir)
    samples_list = sorted(os.listdir("."))
    min_cluster_size = 3
    #bad_regions_file = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/masked_regions/LAN210_v0.10m_masked_all.gff"

    bad_regions_file = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/masked_regions/LAN210_v0.10m_masked_all_not_in_good_genes.gff"
    bad_regions = CollectionGFF(input_file=bad_regions_file,
                                from_file=True)
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

    bad_region_dict = {}
    with open(bad_regions_file) as gff_fd:
        for record in GFF.parse(gff_fd):
            bad_region_dict[record.id] = record
    region = "chrX"
    for sample_set in sample_set_names_list:
        collection = CollectionVCF(from_file=True, vcf_file="%s_good.vcf" % sample_set)
        plt.figure(1, figsize=(15, 8))

        ax0 = plt.subplot(212)
        ax0.set_yscale('log', basey=2)
        ax0.set_xlim(1, 745688)
        #ax0.set_ylim(ymin=0)
        region_reference_dict, maximum = rainfall_plot(collection, region, base_colors=[], facecolor="#D6D6D6",
                      ref_genome=reference, draw_gaps=True, masked_regions=bad_region_dict)
        for reference in region_reference_dict[region]:
            region_reference_dict[region][reference][0] = np.array(region_reference_dict[region][reference][0])
            region_reference_dict[region][reference][1] = np.array(region_reference_dict[region][reference][1])

        ax_list = []
        range_list = [(73400, 74550),
                      (136300, 137500),
                      (159500, 160500),
                      (171800, 173150),
                      (227650, 228550)
                      ]
        for i in range(0, 5):
            ax_list.append(plt.subplot(2, 5, i + 1))
            ax_list[i].set_yscale('log', basey=2)
            ax_list[i].set_xlim(range_list[i][0], range_list[i][1])
            zoom_effect01(ax_list[i], ax0, range_list[i][0], range_list[i][1])
            plt.axhline(y=100, color="#000000")
            plt.axhline(y=1000, color="#000000")
            plt.axhline(y=500, color="purple")
            plt.axhline(y=10, color="#000000")
            if i == 0:
                plt.ylabel("Distance")
            ax_list[i].set_ylim(1, maximum)
            plt.locator_params(axis='x', nbins=2)
            for reference in region_reference_dict[region]:
                in_range = (region_reference_dict[region][reference][0] >= range_list[i][0]) & (region_reference_dict[region][reference][0] <= range_list[i][1])
                x = np.extract(in_range, region_reference_dict[region][reference][0])
                y = np.extract(in_range, region_reference_dict[region][reference][1])
                #print(reference)
                #print(x)
                #print(y)
                plt.plot(x,
                         y,
                         color=reference_colors[reference],
                         marker='.', linestyle='None', label=reference)
        plt.suptitle("Region %s" % region)
        plt.savefig("%s_rainfall_%s.svg" % (sample_set, region))
        plt.savefig("%s_rainfall_%s.png" % (sample_set, region))