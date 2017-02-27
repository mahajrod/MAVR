__author__ = 'mahajrod'

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
os.environ['MPLCONFIGDIR'] = '/tmp/'

import matplotlib.pyplot as plt
plt.ioff()
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
    def draw_alignment(alignment, features, output_prefix, record_style=None, ext_list=["svg", "png"],
                       label_fontsize=13, left_offset=0.2, figure_width=8, id_synonym_dict=None,
                       id_replacement_mode="partial"):
        """
        id_replacement_mode have to be either partial or exact
        """
        from Routines import SequenceRoutines
        sequence_number = len(alignment)
        alignment_length = len(alignment[0].seq)

        figure = plt.figure(figsize=(figure_width, sequence_number))
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
            if (feature.type == "domen") or (feature.type == "domain"):
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
            if id_synonym_dict:
                if id_replacement_mode == "exact":
                    if record.id in id_synonym_dict:
                        record_label = id_synonym_dict[record.id]
                    else:
                        record_label = record.id
                        print("WARNING!!! Synonym for %s was not found" % record.id)
                elif id_replacement_mode == "partial":

                    partial_syn_list = []
                    for partial_syn in id_synonym_dict:
                        if partial_syn in record.id:
                            partial_syn_list.append(partial_syn)

                    if len(partial_syn_list) > 1:
                        print("WARNING!!! More than one possible replacement for %s was found: %s. No replacement then." % (record.id, ",".join(partial_syn_list)))
                        record_label = record.id
                    elif not partial_syn_list:
                        record_label = record.id
                        print("WARNING!!! Synonym for %s was not found" % record.id)
                    else:
                        record_label = id_synonym_dict[partial_syn_list[0]]
                else:
                    raise ValueError("Unknown id replacement mode")

            else:
                record_label = record.id

            subplot.annotate(record_label, xy=(0, gap_y_start), xycoords='data', fontsize=16,
                             xytext=(-15, 1.5 * gap_line_y_shift), textcoords='offset points',
                             ha='right', va='top')

            for gap_coords, gap_len in zip(gap_coords_list, gap_len_list):
                fragment = Rectangle((prev_x, start_y), gap_coords[0] - prev_x - 1, protein_height, fill=False,
                                     edgecolor="black", facecolor="grey")

                subplot.add_patch(fragment)
                prev_x = gap_coords[1]
                plt.plot([gap_coords[0], gap_coords[0] + int(gap_len/2), gap_coords[1] - 1],
                         [gap_y_start, gap_y_jump, gap_y_start], color="black", linewidth=1)

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
                                 xycoords='data', fontsize=label_fontsize,
                                 xytext=(0, 1.5 * gap_line_y_shift), textcoords='offset points', ha='center', va='top')

        plt.xlim(xmax=alignment_length + 10)
        plt.ylim(ymin=0, ymax=start_y + 2 * protein_height)
        #plt.tight_layout()
        plt.subplots_adjust(left=left_offset, right=0.95)#bottom=0.1, right=0.8, top=0.9)
        for extension in ext_list:
            plt.savefig("%s.%s" % (output_prefix, extension))

    def draw_alignment_from_file(self, alignment_file, feature_gff, output_prefix, alignment_style=None,
                                 alignment_format="fasta", ext_list=["svg", "png"], label_fontsize=13,
                                 left_offset=0.2, figure_width=8, id_synonym_dict=None,
                                 id_replacement_mode="partial"):

        alignment = AlignIO.read(alignment_file, format=alignment_format)
        if feature_gff:
            with open(feature_gff, "r") as gff_fd:
                features = list(GFF.parse(gff_fd))[0].features
        else:
            features = []

        self.draw_alignment(alignment, features, output_prefix, ext_list=ext_list, label_fontsize=label_fontsize,
                            left_offset=left_offset, figure_width=figure_width, id_synonym_dict=id_synonym_dict,
                            id_replacement_mode=id_replacement_mode)

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

    def draw_heatmap_and_three_percent_histograms(self, first_histo_values, second_histo_values,
                                                  third_histo_values, output_prefix, figsize=(12, 12),
                                                  extensions=("png", "svg")):
        """
        second_histo_values and third_histo_values are used to build heatmap
        """

        plt.figure(1, figsize=figsize)

        for (index, histo_values, title) in zip([1, 2, 3],
                                                [first_histo_values, second_histo_values, third_histo_values],
                                                ["Total support", "CDS support", "Intron support"]):
            subplot = plt.subplot(2, 2, index)

            self.percent_histogram(histo_values, output_prefix=None, n_bins=20, title=title, xlabel="%",
                                   ylabel="Number of transcripts", label=None, extensions=("png", "svg"),
                                   legend=None, legend_location="best", input_mode="percent", xmax=None,
                                   xmin=None)

        bins = np.linspace(0, 100, 21)

        subplot = plt.subplot(2, 2, 4)
        print bins
        counts, xedges, yedges, image = plt.hist2d(second_histo_values,
                                                   third_histo_values,
                                                   bins=(bins, bins),
                                                   range=[[0, 100], [0, 100]])
        max_counts = int(np.nanmax(counts))

        cmap = plt.get_cmap('jet', max_counts)
        #cmap.set_under('gray')
        mappable = plt.cm.ScalarMappable(cmap=cmap)
        mappable.set_array([])
        mappable.set_clim(0.00001, max_counts)
        #mappable.set_array([])
        #mappable.set_clim(-0.5, ncolors+0.5)
        colorbar = plt.colorbar(mappable)
        plt.xlabel("CDS support")
        plt.ylabel("Intron support")
        plt.title("Transcript support")
        plt.tight_layout()
        for ext in extensions:
            plt.savefig("%s.%s" % (output_prefix, ext))

    @staticmethod
    def draw_precalculated_heatmap(heatmap_array, output_prefix=None, figsize=(5, 5), extensions=("png", "svg")):

        if output_prefix:
            figure = plt.figure(1, figsize=figsize)

        heatmap = plt.imshow(heatmap_array, origin='low', interpolation='none')

        if output_prefix:
            for ext in extensions:
                plt.savefig("%s.%s" % (output_prefix, ext))