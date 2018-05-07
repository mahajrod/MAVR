__author__ = 'mahajrod'
import os
from collections import Iterable, OrderedDict
from math import sqrt

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
import numpy as np

from CustomCollections.GeneralCollections import TwoLvlDict, IdList, SynDict
from Parsers.Abstract import Record, Collection
from Routines import FileRoutines
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

        self.unique_pos_matrix_count, self.unique_pos_matrix_count_percent = self.count_unique_positions_per_sequence()

        #self.unique_insertion_counts, self.unique_gap_counts, self.unique_leading_counts, self.unique_trailing_counts = self.count_unique_positions_per_sequence()

    def __str__(self):
        """
        General stat string

        """

        return "%s\t%i\t%i\t%i\t%f\t%i\t%f\t%i\t%f\t%i\t%f\n" % (self.id, self.number_of_sequences, self.length,
                                                                 max(self.unique_pos_matrix_count[0, ]),
                                                                 max(self.unique_pos_matrix_count_percent[0, ]),
                                                                 max(self.unique_pos_matrix_count[1, ]),
                                                                 max(self.unique_pos_matrix_count_percent[1, ]),
                                                                 max(self.unique_pos_matrix_count[2, ]),
                                                                 max(self.unique_pos_matrix_count_percent[2, ]),
                                                                 max(self.unique_pos_matrix_count[3, ]),
                                                                 max(self.unique_pos_matrix_count_percent[3, ]),)

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
            print "%i sequences in alignment" % number_of_sequences
            print "%i columns in alignment" % alignment_length

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

    def count_unique_positions_per_sequence(self):
        unique_matrix_count = []
        unique_matrix_count_percent = []
        unique_insertion_count_dict = SynDict()
        unique_gap_count_dict = SynDict()

        unique_leading_count_dict = SynDict()
        unique_trailing_count_dict = SynDict()

        for row in range(0, self.number_of_sequences):
            sequence_id = self.alignment[row].id
            unique_gaps = 0
            unique_insertions = 0
            leading_unique_positions = 0
            trailing_unique_positions = 0
            for column in range(0, self.length):
                if self.position_presence_matrix[row, column] == 1:
                    unique_insertions += 1
                elif self.position_presence_matrix[row, column] == -1:
                    unique_gaps += 1

            if (self.position_presence_matrix[row, 0] == 1) or (self.position_presence_matrix[row, 0] == -1):
                for column in range(0, self.length):
                    if self.position_presence_matrix[row, 0] == self.position_presence_matrix[row, column]:
                        leading_unique_positions += self.position_presence_matrix[row, column]
                    else:
                        break

            if (self.position_presence_matrix[row, -1] == 1) or (self.position_presence_matrix[row, -1] == -1):

                for column in range(-1, -self.length):
                    if self.position_presence_matrix[row, -1] == self.position_presence_matrix[row, column]:
                        trailing_unique_positions += self.position_presence_matrix[row, column]
                    else:
                        break

            unique_insertion_count_dict[sequence_id] = unique_insertions
            unique_gap_count_dict[sequence_id] = unique_gaps
            unique_leading_count_dict[sequence_id] = leading_unique_positions
            unique_trailing_count_dict[sequence_id] = trailing_unique_positions

            unique_matrix_count.append([unique_insertions,
                                        unique_gaps,
                                        leading_unique_positions,
                                        trailing_unique_positions])

            unique_matrix_count_percent.append([float(unique_insertions)/float(self.seq_lengths[sequence_id]),
                                                float(unique_gaps)/float(self.seq_lengths[sequence_id]),
                                                float(leading_unique_positions)/float(self.seq_lengths[sequence_id]),
                                                float(trailing_unique_positions)/float(self.seq_lengths[sequence_id])])

        unique_matrix_count = np.array(unique_matrix_count, dtype=int)
        unique_matrix_count_percent = np.array(unique_matrix_count_percent, dtype=float)

        return unique_matrix_count, unique_matrix_count_percent

        #return unique_insertion_count_dict, unique_gap_count_dict, unique_leading_count_dict, unique_trailing_count_dict

    def __str__(self):
        pass


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
    def __init__(self, metadata=None, record_list=None, header=None, input_file=None, filetype="fasta"):
        # metadata should be Metadata class
        # header should be Header class
        # record_list should be list of Record class
        if not(input_file is None):
            self.records = OrderedDict()
            self.read(input_file, filetype=filetype)
        else:
            self.metadata = metadata
            self.records = record_list
            self.header = header

    def read(self, input_file, filetype="fasta"):
        list_of_files = FileRoutines.make_list_of_path_to_files(input_file)
        for filename in list_of_files:
            directory, basename, extension = FileRoutines.split_filename(filename)
            self.records[basename] = MultipleAlignmentStatRecord(basename, alignment=AlignIO.read(input_file, filetype))
        # collectiontype-dependent function
        pass

    def write_general_stats(self, output):
        header = "#aln_id\tseq_number\taln_length\t" \
                 "max_uniq_insertions\tmax_uniq_insertions,%\t" \
                 "max_unique_gaps\tmax_unique_gaps,%\t" \
                 "max_leaiding_unique_pos\tmax_leaiding_unique_pos,%\t" \
                 "max_trailing_unique_pos\tmax_trailing_unique_pos,%\n"

        with open(output, "w") as out_fd:
            out_fd.write(header)
            for record_id in self.records:
                out_fd.write(str(self.records[record_id]) + "\n")
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
