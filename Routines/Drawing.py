__author__ = 'mahajrod'

import os
from collections import OrderedDict

import numpy as np
import matplotlib
matplotlib.use('Agg')
os.environ['MPLCONFIGDIR'] = '/tmp/'

import matplotlib.pyplot as plt
plt.ioff()
from matplotlib.patches import Rectangle, Circle
from matplotlib import text
from BCBio import GFF

from Bio import AlignIO

from CustomCollections.GeneralCollections import SynDict
from Routines.Matplotlib import MatplotlibRoutines
from Routines.Sequence import SequenceRoutines
from Pictures.Features import RectangularProtein
from Parsers.DESeq2 import CollectionPWC


class DrawingRoutines(MatplotlibRoutines, SequenceRoutines):
    def __init__(self):
        MatplotlibRoutines.__init__(self)

    def draw_chromosomes_with_features_simple(self, chromosomes_gff, genes_gff, output_prefix, figsize=(10, 10),
                                              sense_feature_color="green", antisense_feature_color="red",
                                              chromosome_color="black", label_fontsize=15,
                                              chromosome_name_field=None,
                                              ext_list=("png",), dpi=None,
                                              alias_fields_in_gff=("ID", "Name", "gene", "Alias"),
                                              id_field_in_gff="ID", deseq2_pwc_file=None, upregulated_color="green",
                                              downregulated_color="red", absent_expression_data_color="black",
                                              coloring_mode="strand"):
        if deseq2_pwc_file and (coloring_mode == "expression"):
            deseq2_pwc_collection = CollectionPWC(from_file=True, pwc_file=deseq2_pwc_file)

        figure = plt.figure(figsize=figsize, dpi=dpi)
        """
        if dpi:
            figure = plt.figure(figsize=figsize, dpi=dpi)
        else:
            figure = plt.figure(figsize=figsize)
        """
        subplot = plt.subplot(1, 1, 1)

        subplot.get_yaxis().set_visible(False)
        #subplot.get_xaxis().set_visible(False)
        #axes.xaxis.set_major_formatter(x_formatter)

        #subplot.spines['bottom'].set_color('none')
        subplot.spines['right'].set_color('none')
        subplot.spines['left'].set_color('none')
        subplot.spines['top'].set_color('none')

        chr_dict = OrderedDict()

        max_chr_len = 0
        chr_number = 0
        with open(chromosomes_gff, "r") as chr_fd:
            for line in chr_fd:
                if line[0] == "#":
                    continue
                line_list = line.split("\t")
                feature_scaffold = line_list[0]
                feature_type = line_list[2]
                feature_start = line_list[3]
                feature_end = line_list[4]
                feature_description_list = line_list[8].split(";")

                chr_name = None
                if chromosome_name_field:

                    for entry in feature_description_list:
                        entry_list = entry.split("=")
                        if entry_list[0] == chromosome_name_field:
                            chr_name = entry_list[1]

                if feature_type == "chromosome" or feature_type == "plasmid" or feature_type == "region":
                    chr_dict[feature_scaffold] = OrderedDict()
                    chr_dict[feature_scaffold]["Name"] = chr_name if chr_name else feature_scaffold
                    chr_dict[feature_scaffold]["Genes"] = OrderedDict()
                    chr_dict[feature_scaffold]["Length"] = int(feature_end)
                    if chr_dict[feature_scaffold]["Length"] > max_chr_len:
                        max_chr_len = chr_dict[feature_scaffold]["Length"]
                        chr_number += 1
                elif feature_type == "centromere":
                    # add check for case when centromere line is first in gff file

                    chr_dict[feature_scaffold]["Centromere"] = [int(feature_start), int(feature_end)]
        print chr_dict
        if genes_gff:
            with open(genes_gff, "r") as gene_fd:
                for line in gene_fd:
                    if line[0] == "#":
                        continue
                    line_list = line.strip().split("\t")
                    feature_scaffold = line_list[0]
                    feature_type = line_list[2]
                    feature_start = line_list[3]
                    feature_end = line_list[4]
                    feature_strand = line_list[6]
                    description_list = line_list[-1].split(";")

                    entry_id = None
                    alias_list = []
                    print "AAAAAAAAAAA"
                    print line
                    for entry in description_list:
                        if entry[:len(id_field_in_gff)+1] == ("%s=" % id_field_in_gff):
                            entry_id = entry[len(id_field_in_gff)+1:]
                            alias_list.append(entry_id)
                        for alias_id in alias_fields_in_gff:
                            if entry[:len(alias_id)+1] == ("%s=" % alias_id):
                                alias_list += entry.split("=")[1].split(",")
                    print entry_id
                    if not entry_id:
                        continue

                    chr_dict[feature_scaffold]["Genes"][entry_id] = [int(feature_start), int(feature_end), feature_strand, alias_list]

        print chr_dict
        centromere_radius = int(max_chr_len/100)
        distance_between_chromosome_and_gene = centromere_radius * 2
        distance_between_chromosomes = centromere_radius * 8
        chromosome_width = int(centromere_radius / 2)
        gene_width = chromosome_width
        chromosome_position = - int(distance_between_chromosomes / 3)

        text_x_offset = -max_chr_len/15
        for chromosome in chr_dict:
            print "Drawing chromosome %s" % chr_dict[chromosome]["Name"]
            chromosome_position += distance_between_chromosomes
            chromosome_fragments_list = []
            if "Centromere" not in chr_dict[chromosome]:
                chromosome_fragments_list.append(Rectangle((1, chromosome_position), chr_dict[chromosome]["Length"],
                                                           chromosome_width, fill=True, edgecolor=chromosome_color,
                                                           facecolor=chromosome_color))
            else:
                # left arm of chromosome
                print "Centromere"
                #print "Left arm"
                #print (1, chromosome_position), chr_dict[chromosome]["Centromere"][0], chromosome_width


                chromosome_fragments_list.append(Rectangle((1, chromosome_position), chr_dict[chromosome]["Centromere"][0],
                                                           chromosome_width, fill=True, edgecolor=chromosome_color,
                                                           facecolor=chromosome_color))
                chromosome_fragments_list.append(Circle((chr_dict[chromosome]["Centromere"][0] + centromere_radius,
                                                         chromosome_position + chromosome_width/2), centromere_radius,
                                                        fill=False, edgecolor=chromosome_color,
                                                        facecolor=chromosome_color))
                chromosome_fragments_list.append(Rectangle((chr_dict[chromosome]["Centromere"][0] + 2 * centromere_radius,
                                                            chromosome_position), chr_dict[chromosome]["Length"] - chr_dict[chromosome]["Centromere"][1],
                                                           chromosome_width, fill=True, edgecolor=chromosome_color,
                                                           facecolor=chromosome_color))
            for patch in chromosome_fragments_list:
                subplot.add_patch(patch)

            subplot.annotate(chr_dict[chromosome]["Name"], xy=(text_x_offset, chromosome_position),
                             xycoords='data', fontsize=label_fontsize,
                             ha='center', va='center')

            if chr_dict[chromosome]["Genes"]:
                for gene in chr_dict[chromosome]["Genes"]:
                    print "Adding feature %s: %i-%i, %s. Aliases: %s" % (gene, chr_dict[chromosome]["Genes"][gene][0],
                                                                         chr_dict[chromosome]["Genes"][gene][1],
                                                                         chr_dict[chromosome]["Genes"][gene][2],
                                                                         ",".join(chr_dict[chromosome]["Genes"][gene][3]))
                    gene_start = chr_dict[chromosome]["Genes"][gene][0]
                    if "Centromere" in chr_dict[chromosome]:
                        if gene_start >= chr_dict[chromosome]["Centromere"][1]:
                            gene_start += centromere_radius * 2

                    if deseq2_pwc_file and (coloring_mode == "expression"):
                        for alias in chr_dict[chromosome]["Genes"][gene][3]:
                            if alias in deseq2_pwc_collection:
                                gene_log2foldchange = deseq2_pwc_collection[alias].log2foldchange
                                break
                        else:
                            gene_log2foldchange = None
                        gene_color = absent_expression_data_color if gene_log2foldchange is None else upregulated_color if gene_log2foldchange > 0 else downregulated_color
                    else:
                        gene_color = sense_feature_color if chr_dict[chromosome]["Genes"][gene][2] == "+" else antisense_feature_color

                    gene_patch = Rectangle((gene_start, chromosome_position + (1 if chr_dict[chromosome]["Genes"][gene][2] == "+" else -1) * distance_between_chromosome_and_gene),
                                           chr_dict[chromosome]["Genes"][gene][1] - chr_dict[chromosome]["Genes"][gene][0] + 1,
                                           gene_width, fill=True, edgecolor=gene_color,
                                           facecolor=gene_color)
                    subplot.add_patch(gene_patch)
        plt.xlim(xmax=max_chr_len, xmin=text_x_offset)
        plt.ylim(ymax=chromosome_position+distance_between_chromosomes)
        plt.subplots_adjust(right=0.95)#bottom=0.1, right=0.8, top=0.9)
        for extension in ext_list:
            plt.savefig("%s.%s" % (output_prefix, extension), dpi=dpi)

    def draw_alignment(self, alignment, features, output_prefix, record_style=None, ext_list=["svg", "png"],
                       label_fontsize=13, left_offset=0.2, figure_width=8, id_synonym_dict=None,
                       id_replacement_mode="partial", domain_style="vlines"):
        """
        id_replacement_mode have to be either partial or exact
        """
        #from Routines import SequenceRoutines
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

        dist_between_proteins = 10
        start_x = 0
        start_y = - dist_between_proteins

        gap_line_y_shift = int(protein_height/2)
        gap_line_y_jump = int(protein_height/2)

        domen_colors = []
        for feature in features:
            if (feature.type == "domen") or (feature.type == "domain"):
                domen_colors.append(subplot._get_lines.color_cycle.next())

        for record in alignment:

            gap_coords_list, gap_len_list = self.find_homopolymers(record.seq, "-", min_size=1,
                                                                   search_type="perfect")
            #print gap_coords_list, gap_len_list

            start_y += protein_height + dist_between_proteins
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

                if gap_coords[0] != 0:

                    fragment = Rectangle((prev_x, start_y), gap_coords[0] - prev_x, protein_height, fill=False,
                                         edgecolor="black", facecolor="grey")
                    #print prev_x
                    #print gap_coords[0] - prev_x

                    subplot.add_patch(fragment)
                prev_x = gap_coords[1]
                #print [gap_coords[0], gap_coords[0] + int(gap_len/2) + 1, gap_coords[1]]
                plt.plot([gap_coords[0], gap_coords[0] + int(gap_len/2) + 1, gap_coords[1]], #plt.plot([gap_coords[0] + 2, gap_coords[0] + int(gap_len/2) + 1, gap_coords[1] - 1],
                         [gap_y_start, gap_y_jump, gap_y_start], color="black", linewidth=1)

            if not gap_coords_list:
                fragment = Rectangle((prev_x, start_y), alignment_length, protein_height, fill=False,
                                     edgecolor="black", facecolor="grey")
                subplot.add_patch(fragment)
            else:
                if gap_coords_list[-1][-1] != alignment_length:
                    fragment = Rectangle((prev_x, start_y), alignment_length - prev_x, protein_height, fill=False,
                                         edgecolor="black", facecolor="grey")
                    #print prev_x, alignment_length - prev_x
                    subplot.add_patch(fragment)
            i = 0
            for feature in features:
                if (feature.type == "domen") or (feature.type == "domain"):
                    print feature.id, feature.location
                    if domain_style == "rectangle":
                        fragment = Rectangle((feature.location.start, start_y), len(feature)-1, protein_height, fill=False,
                                             facecolor="grey", edgecolor=domen_colors[i]) #edgecolor="green",
                        subplot.add_patch(fragment)
                    elif domain_style == "vlines":
                        plt.vlines(feature.location.start, protein_height, start_y + protein_height, colors=domen_colors[i])
                        plt.vlines(feature.location.end - 1, protein_height, start_y + protein_height, colors=domen_colors[i])
                    i += 1

        for feature in features:
            if feature.type == "domen":
                print feature.id, feature.location
                subplot.annotate(feature.id, xy=(feature.location.start + len(feature)/2, gap_y_start + protein_height),
                                 xycoords='data', fontsize=label_fontsize,
                                 xytext=(0, 1.5 * gap_line_y_shift), textcoords='offset points', ha='center', va='top')

        plt.xlim(xmin=0, xmax=alignment_length + 10)
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
                record = list(GFF.parse(gff_fd))[0]
                features = record.features
                record_id = record.id
                gap_coords_list, gap_len_list = self.find_homopolymers(record.seq, "-", min_size=1,
                                                                       search_type="perfect")
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

    def draw_histogram_from_multiple_files(self, filelist, output_prefix,  filelabels=None, nbins=30,
                                           figsize=(5, 5), title=None, xlabel=None, ylabel=None,
                                           extensions=("png", "svg"), separator="\n", xmin=None, xmax=None,
                                           histtype="stepfilled"):

        figure = plt.figure(1, figsize=figsize)
        data = []

        if filelabels:
            labels = filelabels
        else:
            labels = []
            for filename in filelist:
                path, basename, extention = self.split_filename(filename)
                labels.append(basename)

        if (xmin is not None) and (xmax is not None):
            bins = np.linspace(xmin, xmax, nbins + 1)
        else:
            bins = nbins

        for filename, label in zip(filelist, labels):
            filedata = np.fromfile(filename, sep=separator)
            data.append(filedata)

            print("%s\tmin %f\t max %f" % (label, np.min(filedata), np.max(filedata)))

            a = np.histogram(filedata, bins=bins)
            print a

        print xmin
        print xmax
        print bins

        colors = plt.cm.jet(np.linspace(0, 1, len(filelist)))

        plt.hist(data, label=labels, bins=bins, histtype=histtype, color=colors)

        if xlabel:
            plt.xlabel(xlabel)
        if ylabel:
            plt.ylabel(ylabel)
        if title:
            plt.title(title)

        plt.xlim(xmin=xmin, xmax=xmax)
        plt.legend(loc="best")

        for ext in extensions:
            plt.savefig("%s.%s" % (output_prefix, ext))


















