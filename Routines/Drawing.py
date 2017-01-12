__author__ = 'mahajrod'

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
os.environ['MPLCONFIGDIR'] = '/tmp/'

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib import text
from BCBio import GFF

from Bio import AlignIO

from CustomCollections.GeneralCollections import SynDict
from Routines.Matplotlib import MatplotlibRoutines
from Pictures.Features import RectangularProtein


class DrawingRoutines(MatplotlibRoutines):
    def __init__(self):
        MatplotlibRoutines.__init__(self)

    @staticmethod
    def draw_alignment(alignment, features, output_prefix, record_style=None, ext_list=["svg", "png"]):
        from Routines import SequenceRoutines
        sequence_number = len(alignment)
        alignment_length = len(alignment[0].seq)

        figure = plt.figure(figsize=(8, sequence_number))
        subplot = plt.subplot(1, 1, 1)

        subplot.get_yaxis().set_visible(False)
        #subplot.get_xaxis().set_visible(False)
        #axes.xaxis.set_major_formatter(x_formatter)

        #subplot.spines['bottom'].set_color('none')
        subplot.spines['right'].set_color('none')
        subplot.spines['left'].set_color('none')
        subplot.spines['top'].set_color('none')

        protein_height = 10

        dist_bettwen_proteins = 10
        start_x = 0
        start_y = - dist_bettwen_proteins

        gap_line_y_shift = int(protein_height/2)
        gap_line_y_jump = int(protein_height/2)

        domen_colors = []
        for feature in features:
            if feature.type == "domen":
                domen_colors.append(subplot._get_lines.color_cycle.next())

        for record in alignment:
            print record.id
            gap_coords_list, gap_len_list = SequenceRoutines.find_homopolymers(record.seq, "-", min_size=1,
                                                                               search_type="perfect")

            start_y += protein_height + dist_bettwen_proteins
            gap_y_start = gap_line_y_shift + start_y
            gap_y_jump = gap_y_start + gap_line_y_jump
            prev_x = 0
            """
            figure.text(0, start_y, record.id, rotation=0, fontweight="bold", transform=subplot.transAxes, fontsize=9,
                         horizontalalignment='center',
                         verticalalignment='center')
            """

            subplot.annotate(record.id, xy=(0, gap_y_start), xycoords='data', fontsize=16,
                xytext=(-15, 1.5 * gap_line_y_shift), textcoords='offset points',
                ha='right', va='top')

            for gap_coords, gap_len in zip(gap_coords_list, gap_len_list):
                fragment = Rectangle((prev_x, start_y), gap_coords[0] - prev_x - 1, protein_height, fill=False,
                                     edgecolor="black", facecolor="grey")

                subplot.add_patch(fragment)
                prev_x = gap_coords[1]
                plt.plot([gap_coords[0], gap_coords[0] + int(gap_len/2), gap_coords[1] - 1],
                         [gap_y_start, gap_y_jump, gap_y_start], color="black")

            if not gap_coords_list:
                fragment = Rectangle((prev_x, start_y), alignment_length, protein_height, fill=False,
                                     edgecolor="black", facecolor="grey")
                subplot.add_patch(fragment)
            else:
                if gap_coords_list[-1][-1] != alignment_length:
                    fragment = Rectangle((prev_x, start_y), alignment_length - prev_x, protein_height, fill=False,
                                         edgecolor="black", facecolor="grey")
                    subplot.add_patch(fragment)
            i = 0
            for feature in features:
                if feature.type == "domen":
                    print feature.id, feature.location

                    fragment = Rectangle((feature.location.start, start_y), len(feature)-1, protein_height, fill=False,
                                     facecolor="grey", edgecolor=domen_colors[i]) #edgecolor="green",
                    subplot.add_patch(fragment)
                    i += 1

        for feature in features:
            if feature.type == "domen":
                print feature.id, feature.location
                subplot.annotate(feature.id, xy=(feature.location.start + len(feature)/2, gap_y_start + protein_height),
                                 xycoords='data', fontsize=13,
                                 xytext=(0, 1.5 * gap_line_y_shift), textcoords='offset points', ha='center', va='top')

        plt.xlim(xmax=alignment_length + 1)
        plt.ylim(ymin=0, ymax=start_y + 2*protein_height)
        for extension in ext_list:
            plt.savefig("%s.%s" % (output_prefix, extension))

    def draw_alignment_from_file(self, alignment_file, feature_gff, output_prefix, alignment_style=None,
                                 alignment_format="fasta", ext_list=["svg", "png"]):

        alignment = AlignIO.read(alignment_file, format=alignment_format)

        with open(feature_gff, "r") as gff_fd:
            features = list(GFF.parse(gff_fd))[0].features

        self.draw_alignment(alignment, features, output_prefix, ext_list=ext_list)

    @staticmethod
    def draw_length_histogram(sequence_dict, output_prefix, number_of_bins=None, width_of_bins=None,
                              min_length=1, max_length=None, extensions=("png", "svg"),
                              legend_location='best'):
        length_dict = SynDict()

        for record in sequence_dict:
            length_dict[record] = len(sequence_dict[record].seq)

        length_dict.write("%s.len" % output_prefix)

        lengths = length_dict.values()

        max_len = max(lengths)
        min_len = min(lengths)
        median = np.median(lengths)
        mean = np.mean(lengths)

        if max_length is None:
            maximum_length = max_len
        else:
            maximum_length = max_length

        filtered = []

        if (maximum_length < max_len) and (min_length > 1):
            for entry in lengths:
                if min_length <= entry <= maximum_length:
                    filtered.append(entry)
        elif min_length > 1:
            for entry in lengths:
                if min_length <= entry:
                    filtered.append(entry)
        elif maximum_length < max_len:
            for entry in lengths:
                if entry <= maximum_length:
                    filtered.append(entry)
        else:
            filtered = lengths

        plt.figure(1, figsize=(6, 6))
        plt.subplot(1, 1, 1)

        if number_of_bins:
            bins = number_of_bins
        elif width_of_bins:
            bins = np.arange(min_length - 1, maximum_length, width_of_bins, dtype=np.int32)
            bins[0] += 1
            bins = np.append(bins, [maximum_length])
        else:
            bins = 30
        plt.hist(filtered, bins=bins)
        plt.xlim(xmin=min_length, xmax=maximum_length)
        plt.xlabel("Length")
        plt.ylabel("N")
        plt.title("Distribution of sequence lengths")
        plt.legend(("Min: %i\nMax: %i\nMean: %i\nMedian: %i" % (min_len, max_len, mean, median),), loc=legend_location)
        for ext in extensions:
            plt.savefig("%s.%s" % (output_prefix, ext))

        os.remove("temp.idx")
