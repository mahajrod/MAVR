__author__ = 'mahajrod'
import os
import math
from collections import Iterable, OrderedDict
from math import sqrt

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
import numpy as np

from CustomCollections.GeneralCollections import TwoLvlDict, IdList, SynDict
from Parsers.Abstract import Record, Collection
from Routines import FileRoutines, MatplotlibRoutines
from Bio import AlignIO


class MultipleAlignmentStatRecord(Record):

    def __init__(self, aln_id, alignment, description=None, flags=None, gap_symbol="-", verbose=False):
        self.id = aln_id
        self.alignment = alignment
        self.description = description
        self.flags = flags
        self.gap_symbol = gap_symbol
        self.number_of_sequences = len(alignment)
        self.length = alignment.get_alignment_length()

        self.position_presence_matrix = self.get_position_presence_matrix(verbose=verbose)
        self.gap_counts, self.seq_lengths = self.count_gaps()
        self.N_counts = self.count_Ns()

        self.unique_pos_matrix_count, self.unique_pos_matrix_count_percent = self.count_unique_positions_per_sequence()

        #self.unique_insertion_counts, self.unique_gap_counts, self.unique_leading_counts, self.unique_trailing_counts = self.count_unique_positions_per_sequence()

    def __str__(self):
        """
        General stat string

        """
        """
        return "%s\t%i\t%i\t%i\t%i\t%i\t%f\t%i\t%f\t%i\t%f\t%i\t%f\t%i\t%f\n" % (self.id, self.number_of_sequences, self.length,
                                                                         min(self.seq_lengths.values()),
                                                                         max(self.seq_lengths.values()),
                                                                         max(self.unique_pos_matrix_count[:, 0]),
                                                                         max(self.unique_pos_matrix_count_percent[:, 0]),
                                                                         max(self.unique_pos_matrix_count[:, 1]),
                                                                         max(self.unique_pos_matrix_count_percent[:, 1]),
                                                                         max(self.unique_pos_matrix_count[:, 2]) if math.fabs(max(self.unique_pos_matrix_count[:, 2])) > math.fabs(min(self.unique_pos_matrix_count[:, 2])) else min(self.unique_pos_matrix_count[:, 2]),
                                                                         max(self.unique_pos_matrix_count_percent[:, 2]) if math.fabs(max(self.unique_pos_matrix_count_percent[:, 2])) > math.fabs(min(self.unique_pos_matrix_count_percent[:, 2])) else min(self.unique_pos_matrix_count_percent[:, 2]),
                                                                         max(self.unique_pos_matrix_count[:, 3]) if math.fabs(max(self.unique_pos_matrix_count[:, 3])) > math.fabs(min(self.unique_pos_matrix_count[:, 3])) else min(self.unique_pos_matrix_count[:, 3]),
                                                                         max(self.unique_pos_matrix_count_percent[:, 3]) if math.fabs(max(self.unique_pos_matrix_count_percent[:, 3])) > math.fabs(min(self.unique_pos_matrix_count_percent[:, 3])) else min(self.unique_pos_matrix_count_percent[:, 3]),
                                                                         max(self.unique_pos_matrix_count[:, 4]),
                                                                         max(self.unique_pos_matrix_count_percent[:, 4]))
        """
        return "%s\t%s\n" % (self.id, "\t".join(map(str, self.get_general_stats())))

    def count_gaps(self):

        gaps_dict = SynDict()
        seq_length_dict = SynDict()

        for row in range(0, self.number_of_sequences):
            sequence_id = self.alignment[row].id
            gaps_dict[sequence_id] = str(self.alignment[row].seq).count(self.gap_symbol)
            seq_length_dict[sequence_id] = self.length - gaps_dict[sequence_id]

        return gaps_dict, seq_length_dict

    def get_position_presence_matrix(self, verbose=True):

        number_of_sequences = len(self.alignment)
        alignment_length = len(self.alignment[0])
        # converting alignment to numpy letter array stored by columns!
        alignment_array = np.array([list(rec) for rec in self.alignment], np.character, order="F")

        if verbose:
            print("%i sequences in alignment" % number_of_sequences)
            print("%i columns in alignment" % alignment_length)

        position_presence_array = np.array([[0 for letter in rec.seq] for rec in self.alignment], int, order="F")

        for column in range(0, alignment_length):
            column_list = alignment_array[:, column]
            number_of_gaps = 0
            for element in column_list:
                if element == self.gap_symbol:
                    number_of_gaps += 1
            for row in range(0, number_of_sequences):
                position_presence_array[row, column] = -number_of_gaps if alignment_array[row, column] == self.gap_symbol else number_of_sequences - number_of_gaps

        return position_presence_array

    def count_Ns(self):
        Ns_dict = SynDict()
        for row in range(0, self.number_of_sequences):
            sequence_id = self.alignment[row].id
            Ns_dict[sequence_id] = self.alignment[row].seq.count("N")
        return Ns_dict

    def get_general_stats(self):
        return [self.number_of_sequences, self.length,
                min(self.seq_lengths.values()),
                max(self.seq_lengths.values()),
                max(self.unique_pos_matrix_count[:, 0]),
                max(self.unique_pos_matrix_count_percent[:, 0]),
                max(self.unique_pos_matrix_count[:, 1]),
                max(self.unique_pos_matrix_count_percent[:, 1]),
                max(self.unique_pos_matrix_count[:, 2]) if math.fabs(max(self.unique_pos_matrix_count[:, 2])) > math.fabs(min(self.unique_pos_matrix_count[:, 2])) else min(self.unique_pos_matrix_count[:, 2]),
                max(self.unique_pos_matrix_count_percent[:, 2]) if math.fabs(max(self.unique_pos_matrix_count_percent[:, 2])) > math.fabs(min(self.unique_pos_matrix_count_percent[:, 2])) else min(self.unique_pos_matrix_count_percent[:, 2]),
                max(self.unique_pos_matrix_count[:, 3]) if math.fabs(max(self.unique_pos_matrix_count[:, 3])) > math.fabs(min(self.unique_pos_matrix_count[:, 3])) else min(self.unique_pos_matrix_count[:, 3]),
                max(self.unique_pos_matrix_count_percent[:, 3]) if math.fabs(max(self.unique_pos_matrix_count_percent[:, 3])) > math.fabs(min(self.unique_pos_matrix_count_percent[:, 3])) else min(self.unique_pos_matrix_count_percent[:, 3]),
                max(self.unique_pos_matrix_count[:, 4]),
                max(self.unique_pos_matrix_count_percent[:, 4])]

    def count_unique_positions_per_sequence(self):
        unique_matrix_count = []
        unique_matrix_count_percent = []
        #unique_insertion_count_dict = SynDict()
        #unique_gap_count_dict = SynDict()

        #unique_leading_count_dict = SynDict()
        #unique_trailing_count_dict = SynDict()

        for row in range(0, self.number_of_sequences):
            sequence_id = self.alignment[row].id
            unique_gaps = 0
            unique_insertions = 0
            leading_unique_positions = 0
            trailing_unique_positions = 0
            unique_positions = 0
            for column in range(0, self.length):
                if self.position_presence_matrix[row, column] == 1:
                    unique_insertions += 1
                    unique_positions += 1
                elif self.position_presence_matrix[row, column] == -1:
                    unique_gaps += 1
                    unique_positions += 1

            if (self.position_presence_matrix[row, 0] == 1) or (self.position_presence_matrix[row, 0] == -1):
                for column in range(0, self.length):
                    if self.position_presence_matrix[row, 0] == self.position_presence_matrix[row, column]:
                        leading_unique_positions += self.position_presence_matrix[row, column]
                    else:
                        break

            if (self.position_presence_matrix[row, -1] == 1) or (self.position_presence_matrix[row, -1] == -1):

                for column in range(-1, -self.length, -1):
                    if self.position_presence_matrix[row, -1] == self.position_presence_matrix[row, column]:
                        trailing_unique_positions += self.position_presence_matrix[row, column]
                    else:
                        break

            #unique_insertion_count_dict[sequence_id] = unique_insertions
            #unique_gap_count_dict[sequence_id] = unique_gaps
            #unique_leading_count_dict[sequence_id] = leading_unique_positions
            #unique_trailing_count_dict[sequence_id] = trailing_unique_positions
            #print unique_trailing_count_dict[sequence_id]

            unique_matrix_count.append([unique_insertions,
                                        unique_gaps,
                                        leading_unique_positions,
                                        trailing_unique_positions,
                                        unique_positions])

            unique_matrix_count_percent.append([float(unique_insertions)/float(self.seq_lengths[sequence_id]),
                                                float(unique_gaps)/float(self.seq_lengths[sequence_id]),
                                                float(leading_unique_positions)/float(self.seq_lengths[sequence_id]),
                                                float(trailing_unique_positions)/float(self.seq_lengths[sequence_id]),
                                                float(unique_positions)/float(self.seq_lengths[sequence_id])])

        unique_matrix_count = np.array(unique_matrix_count, dtype=int)
        unique_matrix_count_percent = np.array(unique_matrix_count_percent, dtype=float)

        return unique_matrix_count, unique_matrix_count_percent

        #return unique_insertion_count_dict, unique_gap_count_dict, unique_leading_count_dict, unique_trailing_count_dict


class Metadata():
    def __init__(self, metadata=[]):
        self.metadata = metadata

    def add_metadata(self):
        pass

    def __str__(self):
        pass


class Header():
    def __init__(self):
        pass

    def __str__(self):
        pass


class MultipleAlignmentStatCollection(Collection):
    # TODO: develop this class to minimize new code in various collections
    def __init__(self, metadata=None, record_list=None, header=None, input_file=None, filetype="fasta", verbose=False,
                 general_stats_input=False):
        # metadata should be Metadata class
        # header should be Header class
        # record_list should be list of Record class
        self.general_stats_table = None
        self.general_stats_header = "#aln_id\tseq_number\taln_length\t" \
                                    "min_seq_len\tmax_seq_len\t" \
                                    "max_uniq_insertions\tmax_uniq_insertions,fraction\t" \
                                    "max_unique_gaps\tmax_unique_gaps,fraction\t" \
                                    "max_leading_unique_pos\tmax_leading_unique_pos,fraction\t" \
                                    "max_trailing_unique_pos\tmax_trailing_unique_pos,fraction\t" \
                                    "max_unique_pos\tmax_unique_pos,fraction\n"
        if not(input_file is None):
            self.records = OrderedDict()
            self.record_id_list = []
            if general_stats_input:
                self.read_general_stats_file(input_file)

            else:
                self.read(input_file, filetype=filetype, verbose=verbose)
                self.general_stats_table = self.get_general_stats()
        else:
            self.metadata = metadata
            self.records = record_list
            self.header = header

    def read(self, input_file, filetype="fasta", verbose=False):
        list_of_files = FileRoutines.make_list_of_path_to_files(input_file)
        for filename in list_of_files:
            if verbose:
                print("Parsing %s ..." % filename)
            directory, basename, extension = FileRoutines.split_filename(filename)
            try:
                self.records[basename] = MultipleAlignmentStatRecord(basename, alignment=AlignIO.read(filename, filetype))
                self.record_id_list.append(basename)
            except:
                raise ValueError("ERROR: Issues while parsing or calculating stats for %s!!!" % filename)
        # collectiontype-dependent function
        pass

    def read_general_stats_file(self, general_stats_file):
        self.record_id_list = []
        general_stats_list = []
        with open(general_stats_file, "r") as in_fd:
            in_fd.readline()
            for line in in_fd:
                line_list = line.strip("\n").split("\t")
                self.record_id_list.append(line_list[0])
                for i in [1, 2, 3, 4, 5, 7, 9, 11, 13]:
                    line_list[i] = int(line_list[i])
                for i in [6, 8, 10, 12, 14]:
                    line_list[i] = float(line_list[i])
                general_stats_list.append(line_list[1:])
        self.general_stats_table = np.array(general_stats_list)

    def get_general_stats(self):
        general_stats_list = []
        for record_id in self.records:
            #print self.records
            general_stats_list.append(self.records[record_id].get_general_stats())
        return np.array(general_stats_list)

    def write_general_stats(self, output):
        with open(output, "w") as out_fd:
            out_fd.write(self.general_stats_header)
            for record_id in self.records:
                #print self.records
                out_fd.write(str(self.records[record_id]))

    def write_stats(self, output_prefix):
        Ns_dict = TwoLvlDict()
        gaps_dict = TwoLvlDict()
        for record_id in self.records:
            Ns_dict[self.records[record_id].id] = self.records[record_id].N_counts
            gaps_dict[self.records[record_id].id] = self.records[record_id].gap_counts

        Ns_dict.write(out_filename="%s.N_counts" % output_prefix)
        gaps_dict.write(out_filename="%s.gaps_counts" % output_prefix)

    def draw_general_stats_distributions(self, output_prefix, figsize=(12, 6), extensions=("png", "svg"), dpi=300,
                                         logscale_heatmaps=True):

        nrows = 2
        ncols = 4

        #figure = plt.figure(figsize=(12, 6), dpi=dpi)
        #ax_array = figure.subplots(nrows=nrows, ncols=ncols, squeeze=False)
        percent_histogram_bin_number = 20

        min_seq_number_in_alignment = min(self.general_stats_table[:, 0])
        max_seq_number_in_alignment = max(self.general_stats_table[:, 0])
        seq_bin_number = max_seq_number_in_alignment - min_seq_number_in_alignment

        figure, ax_array = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize, squeeze=False, dpi=dpi)

        # Histogram: distribution of sequence number in alignment
        print("Histogram: distribution of sequence number in alignment")

        MatplotlibRoutines.draw_histogram(self.general_stats_table[:, 0],
                                          output_prefix="%s.seq_number_distibution" % output_prefix,
                                          number_of_bins=None, width_of_bins=1,
                                          max_threshold=None, min_threshold=None, xlabel=None, ylabel="N of alignments",
                                          title="Distribution of\nsequence numbers",
                                          extensions=("png",), ylogbase=None, subplot=ax_array[0, 0], suptitle=None,
                                          close_figure=False, save_histovalues_only=True)
        # Histogram(logscaled): distribution of sequence number in alignment

        print("Histogram(logscaled): distribution of sequence number in alignment")
        MatplotlibRoutines.draw_histogram(self.general_stats_table[:, 0], output_prefix=None,
                                          number_of_bins=None, width_of_bins=1,
                                          max_threshold=None, min_threshold=None,
                                          xlabel="N of sequences\nin alignment", ylabel="N of alignments",
                                          title=None,
                                          extensions=("png",), ylogbase=10, subplot=ax_array[1, 0], suptitle=None,
                                          close_figure=False)

        #print self.general_stats_table[:, 4].astype(float) / self.general_stats_table[:, 2].astype(float) * 100

        # Heatmap: x: max_seq_len/aln_len, y: min_seq_len/aln_len
        print("Heatmap: x: max_seq_len/aln_len, y: min_seq_len/aln_len")
        MatplotlibRoutines.draw_percent_heatmap(self.general_stats_table[:, 3].astype(float) / self.general_stats_table[:, 1].astype(float) * 100,
                                                self.general_stats_table[:, 2].astype(float) / self.general_stats_table[:, 1].astype(float) * 100,
                                                output_prefix="%s.min_max_seq_len" % output_prefix, xlabel="Max seq len, % of aln",
                                                ylabel="Min seq len, % of aln", title=None,
                                                figsize=(8, 8), minimum_counts_to_show=1,
                                                extensions=("png", "svg"), show_colorbar=True,
                                                bin_number=percent_histogram_bin_number, bin_width=None, bin_array=None,
                                                type="percent", add_max_value=True,
                                                subplot=ax_array[0, 1],
                                                header="#left_xedge\tleft_yedge\tvalue",
                                                save_histovalues_only=True,
                                                logscaled=logscale_heatmaps)

        # Heatmap; x: seq number, y: max_unique_pos_in_seq/aln_len
        print("Heatmap: x: seq number, y: max_unique_pos_in_seq/aln_len")
        MatplotlibRoutines.draw_heatmap(self.general_stats_table[:, 0], 100 * self.general_stats_table[:, 13],
                                        output_prefix="%s.max_unique_pos_seq_number" % output_prefix,
                                        xlabel="N of sequences\nin alignment",
                                        ylabel="Max uniq positions, % of seq", title=None,
                                        figsize=figsize, minimum_counts_to_show=1,
                                        extensions=extensions, show_colorbar=True,
                                        bin_number=(seq_bin_number, percent_histogram_bin_number),
                                        bin_width=None, bin_array=None, min_x_value=min_seq_number_in_alignment,
                                        max_x_value=max_seq_number_in_alignment, min_y_value=0, max_y_value=100,
                                        add_max_value=True, subplot=ax_array[1, 1],
                                        save_histovalues_only=True, header="#left_xedge\tleft_yedge\tcounts",
                                        logscaled=logscale_heatmaps)

        # Heatmap; x: seq number, y: max_unique_insertions_in_seq/aln_len
        print("Heatmap: x: seq number, y: max_unique_insertions_in_seq/aln_len")
        MatplotlibRoutines.draw_heatmap(self.general_stats_table[:, 0], 100 * self.general_stats_table[:, 5],
                                        output_prefix="%s.max_unique_insertions_seq_number" % output_prefix,
                                        xlabel=None, #"N of sequences\nin alignment",
                                        ylabel="Max uniq insertions, % of seq", title=None,
                                        figsize=figsize, minimum_counts_to_show=1,
                                        extensions=extensions, show_colorbar=True,
                                        bin_number=(seq_bin_number, percent_histogram_bin_number),
                                        bin_width=None, bin_array=None, min_x_value=min_seq_number_in_alignment,
                                        max_x_value=max_seq_number_in_alignment, min_y_value=0, max_y_value=100,
                                        add_max_value=True, subplot=ax_array[0, 2],
                                        save_histovalues_only=True, header="#left_xedge\tleft_yedge\tcounts",
                                        logscaled=logscale_heatmaps)

        # Heatmap: x: seq number, y: max_unique_gaps_in_seq/aln_len
        print("Heatmap: x: seq number, y: max_unique_gaps_in_seq/aln_len")
        MatplotlibRoutines.draw_heatmap(self.general_stats_table[:, 0], 100 * self.general_stats_table[:, 7],
                                        output_prefix="%s.max_unique_gaps_seq_number" % output_prefix,
                                        xlabel="N of sequences\nin alignment",
                                        ylabel="Max uniq gaps, % of seq", title=None,
                                        figsize=figsize, minimum_counts_to_show=1,
                                        extensions=extensions, show_colorbar=True,
                                        bin_number=(seq_bin_number, percent_histogram_bin_number),
                                        bin_width=None, bin_array=None, min_x_value=min_seq_number_in_alignment,
                                        max_x_value=max_seq_number_in_alignment, min_y_value=0, max_y_value=100,
                                        add_max_value=True, subplot=ax_array[1, 2],
                                        save_histovalues_only=True, header="#left_xedge\tleft_yedge\tcounts",
                                        logscaled=logscale_heatmaps)


        # Heatmap; x: seq number, y: max_unique_leading_pos_in_seq/aln_len
        print("Heatmap: x: seq number, y: max_unique_leading_pos_in_seq/aln_len")
        MatplotlibRoutines.draw_heatmap(self.general_stats_table[:, 0], 100 * self.general_stats_table[:, 9],
                                        output_prefix="%s.max_unique_leading_pos_seq_number" % output_prefix,
                                        xlabel=None, #"N of sequences\nin alignment",
                                        ylabel="Max uniq leading pos, % of seq", title=None,
                                        figsize=figsize, minimum_counts_to_show=1,
                                        extensions=extensions, show_colorbar=True,
                                        bin_number=(seq_bin_number, percent_histogram_bin_number),
                                        bin_width=None, bin_array=None, min_x_value=min_seq_number_in_alignment,
                                        max_x_value=max_seq_number_in_alignment, min_y_value=0, max_y_value=100,
                                        add_max_value=True, subplot=ax_array[0, 3],
                                        save_histovalues_only=True, header="#left_xedge\tleft_yedge\tcounts",
                                        logscaled=logscale_heatmaps)

        # Heatmap: x: seq number, y: max_unique_trailing_pos_in_seq/aln_len
        print("Heatmap: x: seq number, y: max_unique_trailing_pos_in_seq/aln_len")
        MatplotlibRoutines.draw_heatmap(self.general_stats_table[:, 0], 100 * self.general_stats_table[:, 11],
                                        output_prefix="%s.max_unique_trailing_pos_seq_number" % output_prefix,
                                        xlabel="N of sequences\nin alignment",
                                        ylabel="Max uniq trailing pos, % of seq", title=None,
                                        figsize=figsize, minimum_counts_to_show=1,
                                        extensions=extensions, show_colorbar=True,
                                        bin_number=(seq_bin_number, percent_histogram_bin_number),
                                        bin_width=None, bin_array=None, min_x_value=min_seq_number_in_alignment,
                                        max_x_value=max_seq_number_in_alignment, min_y_value=0, max_y_value=100,
                                        add_max_value=True, subplot=ax_array[1, 3],
                                        save_histovalues_only=True, header="#left_xedge\tleft_yedge\tcounts",
                                        logscaled=logscale_heatmaps)

        plt.tight_layout()
        #aaa = self.general_stats_table[:, 3].astype(float) / self.general_stats_table[:, 1].astype(float) * 100
        #for i in range(0, len(aaa)):
        #    print "%i\t%i\t%i\t%f" % (i, self.general_stats_table[:, 3][i], self.general_stats_table[:, 1][i], aaa[i])

        for ext in extensions:
            plt.savefig("%s.%s" % (output_prefix, ext))

    """
    def __str__(self):
        collection_string = ""
        if self.metadata:
            collection_string += str(self.metadata)
        if self.header:
            collection_string += "\n" + str(self.header)
        if self.records:
            collection_string += "\n" + "\n".join([str(record) for record in self.records])

        return collection_string

    def __iter__(self):
        for record in self.records:
            yield record

    def __getitem__(self, item):
        #string-like access
        return self.records[item]

    def __len__(self):
        return len(self.records)

    def pop(self, index=None):
        # i am a fairy idiot - forgot that 0 is equal to False in python.
        # What to do? Everywhere make direct check if variable is None
        return self.records.pop() if index is None else self.records.pop(index)

    def filter_records_by_parameter_value(self, parameter, minimum=None, maximum=None):
        # TODO: check
        filtered_records = []
        filtered_out_records = []
        if (minimum is not None) and (maximum is not None):
            comparison = lambda par, min_value, max_value: min_value <= par <= max_value
        elif (minimum is None) and (maximum is None):
            raise ValueError("Both minimum and maximum thresholds were not set")
        elif minimum is None:
            comparison = lambda par, min_value, max_value: par <= max_value
        else:
            comparison = lambda par, min_value, max_value: min_value <= par

        for record in self.records:
            if comparison(getattr(record, parameter), minimum, maximum):
                filtered_records.append(record)
            else:
                filtered_out_records.append(record)

        return filtered_records, filtered_out_records

    def filter_records_by_expression(self, expression):
        filtered_records = []
        filtered_out_records = []
        for record in self.records:
            if eval(expression):
                filtered_records.append(record)
            else:
                filtered_out_records.append(record)

        return filtered_records, filtered_out_records

    def filter_records_by_function(self, function):
        filtered_records = []
        filtered_out_records = []
        for record in self.records:
            #print("a\na\na\na\na\na\na\na")
            #print(record.description["Power"])
            if function(record):
                filtered_records.append(record)
            else:
                filtered_out_records.append(record)

        return filtered_records, filtered_out_records

    def check_location(self, bad_region_collection_gff, expression="bad_region.start <= self.pos <= bad_region.end"):
        for record in self:
            record.check_location(bad_region_collection_gff, expression=expression)
    """
    """
    def filter_by_expression(self, expression):
        self_type = self.__class__.__name__
        filtered_records, filtered_out_records = self.filter_records_by_expression(expression)
        #TODO: think how to rewrite following damned string!!!!!!!!!!
        exec("from Parsers.%s import %s" % (self_type[10:], self_type))
        return eval(self_type)(metadata=self.metadata, record_list=filtered_records,
                               header_list=self.header_list, from_file=False), \
               eval(self_type)(metadata=self.metadata, record_list=filtered_out_records,
                               header_list=self.header_list, from_file=False)
    """
    """
    def split_records_by_flags(self, flag_set, mode="all"):
        # possible modes:
        # all - record to be counted as 'with flag' must have all flags from flags_list
        # one - record to be counted as 'with flag' must have at least one flag from flags_list
        flags = set(flag_set)
        with_flags_records = []
        without_flags_records = []
        if mode == "all":
            for record in self:
                if flags & record.flags == flags:
                    with_flags_records.append(record)
                else:
                    without_flags_records.append(record)
        elif mode == "one":
            for record in self:
                if flags & record.flags:
                    with_flags_records.append(record)
                else:
                    without_flags_records.append(record)

        return with_flags_records, without_flags_records

    def write(self, output_file):
        with open(output_file, "w") as out_fd:
            if self.metadata:
                out_fd.write(str(self.metadata) + "\n")
            if self.header:
                out_fd.write(str(self.header) + "\n")
            for record in self.records:
                out_fd.write(str(record) + "\n")

    def get_location(self, record_dict, key="Loc", use_synonym=False, synonym_dict=None):
        for record in self.records:
            record.get_location(record_dict, key=key, use_synonym=use_synonym, synonym_dict=synonym_dict)

    def _split_regions(self):
        splited_dict = OrderedDict({})
        for record in self.records:
            if record.chrom not in splited_dict:
                splited_dict[record.chrom] = [record]
            else:
                splited_dict[record.chrom].append(record)
        return splited_dict
    """
