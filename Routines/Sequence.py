__author__ = 'mahajrod'
import os
import re
import sys
from copy import deepcopy
from collections import OrderedDict

import numpy as np

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from CustomCollections.GeneralCollections import TwoLvlDict, SynDict, IdList
from Routines.Functions import output_dict


class SequenceRoutines():

    def __init__(self):
        pass

    @staticmethod
    def get_lengths(record_dict, out_file=None, close_after_if_file_object=False):
        lengths_dict = SynDict()
        for record_id in record_dict:
            lengths_dict[record_id] = len(record_dict[record_id])

        if out_file:
            lengths_dict.write(out_file, header=False, separator="\t", splited_values=False, values_separator=",",
                               close_after_if_file_object=close_after_if_file_object)

        return lengths_dict

    @staticmethod
    def get_lengths_from_seq_file(input_file_list, format="fasta", out_file=None, close_after_if_file_object=False):
        record_dict = SeqIO.index_db("tmp.idx", input_file_list, format=format)
        lengths_dict = SynDict()
        for record_id in record_dict:
            lengths_dict[record_id] = len(record_dict[record_id])

        if out_file:
            lengths_dict.write(out_file, header=False, separator="\t", splited_values=False, values_separator=",",
                               close_after_if_file_object=close_after_if_file_object)

        os.remove("tmp.idx")
        return lengths_dict

    @staticmethod
    def record_by_id_generator(record_dict, id_list, verbose=False):
        for record_id in id_list:
            if record_id in record_dict:
                yield record_dict[record_id]
            else:
                if verbose:
                    sys.stderr.write("Not found: %s\n" % record_id)

    def extract_sequence_by_ids(self, sequence_file, id_file, output_file, format="fasta", verbose=False):
        tmp_index_file = "tmp.idx"
        id_list = IdList()
        id_list.read(id_file)
        if verbose:
            print("Parsing %s..." % (sequence_file if isinstance(id_file, str) else ",".join(id_file)))

        sequence_dict = SeqIO.index_db(tmp_index_file, sequence_file, format=format)
        SeqIO.write(self.record_by_id_generator(sequence_dict, id_list, verbose=verbose),
                    output_file, format=format)
        os.remove(tmp_index_file)

    @staticmethod
    def find_gaps(record_dict):
        gap_reg_exp = re.compile("N+", re.IGNORECASE)
        gaps_dict = {}
        for region in record_dict:
            gaps_dict[region] = SeqRecord(seq=record_dict[region].seq,
                                          id=record_dict[region].id,
                                          description=record_dict[region].description)
            gaps = gap_reg_exp.finditer(str(record_dict[region].seq))  # iterator with
            for match in gaps:
                gaps_dict[region].features.append(SeqFeature(FeatureLocation(match.start(), match.end()),
                                                  type="gap", strand=None))
        return gaps_dict

    @staticmethod
    def record_by_expression_generator(record_dict, expression=None, id_file="passed_records.ids"):
        """
        :param record_dict: dictionary containing Biopython SeqRecords as values
        :param expression: function to apply to all records in record_dict. If it returns True record will be yielded
        :param id_file: file to write ids of records passed expression
        :return: None
        """
        id_fd = open(id_file, "w")
        if expression is None:
            for record_id in record_dict:
                id_fd.write(record_id + "\n")
                yield record_dict[record_id]
        else:
            for record_id in record_dict:
                if expression(record_dict[record_id]):
                    id_fd.write(record_id + "\n")
                    yield record_dict[record_id]
        id_fd.close()

def get_lengths(record_dict, out_file="lengths.t", write=False, write_header=True):
    lengths_dict = OrderedDict({})
    for record_id in record_dict:
        lengths_dict[record_id] = len(record_dict[record_id])

    if write:
        output_dict(lengths_dict, out_file=out_file, write=write,
                    header_tuple=("record", "length") if write_header else None)

    return lengths_dict


def get_feature_lengths(record_dict):
    lengths_dict = OrderedDict({})

    for record_id in record_dict:
        lengths_dict[record_id] = {}
        for feature in record_dict[record_id].features:
            if feature.type not in lengths_dict[record_id]:
                lengths_dict[record_id][feature.type] = []
            lengths_dict[record_id][feature.type].append(len(feature))
            if feature.type == "gene":
                tmp_dict = {}
                for sub_feature in feature.sub_features:
                    if sub_feature.type not in tmp_dict:
                        tmp_dict[sub_feature.type] = 0
                    if sub_feature.type not in lengths_dict[record_id]:
                        lengths_dict[record_id][sub_feature.type] = []
                    tmp_dict[sub_feature.type] += len(sub_feature)

                for sub_feature_type in tmp_dict:
                    lengths_dict[record_id][sub_feature_type].append(tmp_dict[sub_feature_type])
    record_ids = lengths_dict.keys()
    lengths_dict["all"] = {}
    for record_id in record_ids:
        for feature_type in lengths_dict[record_id]:
            if feature_type not in lengths_dict["all"]:
                lengths_dict["all"][feature_type] = []
            lengths_dict["all"][feature_type] += lengths_dict[record_id][feature_type]

    #for record_id in lengths_dict:
    #    for feature_type in lengths_dict[record_id]:
    #        lengths_dict[record_id][feature_type] = np.array(lengths_dict[record_id][feature_type], dtype=int)
    return lengths_dict


def feature_lengths_collapse_records(lengths_dict,
                                     synonym_dict=None):
    tmp_dict = synonym_dict if synonym_dict is not None else {}
    collapsed_dict = {}
    for record_id in lengths_dict:
        if record_id == "all":
            continue
        for feature_type in lengths_dict[record_id]:
            tmp_type = tmp_dict[feature_type] if feature_type in tmp_dict else feature_type
            if tmp_type not in collapsed_dict:
                collapsed_dict[tmp_type] = lengths_dict[record_id][feature_type]
            else:
                collapsed_dict[tmp_type] += lengths_dict[record_id][feature_type]

    return collapsed_dict


def get_total_feature_lengths(lengths_dict, out_filename=None):
    total_lengths = TwoLvlDict({})
    for record_id in lengths_dict:
        total_lengths[record_id] = {}
        for feature_type in lengths_dict[record_id]:
            total_lengths[record_id][feature_type] = sum(lengths_dict[record_id][feature_type])
    if out_filename is not None:
        total_lengths.write(out_filename, absent_symbol="0")

    return total_lengths


# --------------------Search--------------------------
#              !!!!IMPORTANT!!!!!!
# in this section python notation for coordinates inside sequence is used
def find_homopolymer_end(seq, nucleotide, seq_length, start, search_type="perfect",
                       max_single_insert_size=1, max_total_insert_length=None, max_number_of_insertions=2):
    shift_list = []
    number_of_insertions = 0
    total_insertion_length = 0
    insertion_length = 0
    i = start + 1
    if search_type == 'perfect':
        while i < seq_length:
            if seq[i] != nucleotide:
                return i, None
            i += 1
        return seq_length, None
    else:
        while i < seq_length:
            if seq[i] != nucleotide:
                if seq[i-1] == nucleotide:
                    shift_list.append(i)
                    insertion_length = 1
                else:
                    insertion_length += 1
                number_of_insertions += 1
                total_insertion_length += 1
                if number_of_insertions > max_number_of_insertions or insertion_length > max_single_insert_size:
                    end = shift_list[-1]
                    break
                if max_total_insert_length:
                    if total_insertion_length > max_total_insert_length:
                        end = shift_list[-1]
                        break
            i += 1
        else:
            if seq[-1] == nucleotide:
                end = seq_length
            else:
                end = shift_list[0]

        new_start = shift_list[0] if shift_list else end
        return end, new_start


def find_homopolymers(seq, nucleotide, min_size=5, search_type="perfect",
                      max_single_insert_size=1, max_total_insert_length=None, max_number_of_insertions=2):
    # search types:
    #   perfect - search only for perfect homopolymers, all options other than min_size are ignored
    #   non_perfect - search for non_perfect homopolymers with max_single_insert_size, max_total_insert_length
    #                 and max_number_of_insertions
    seq_length = len(seq)
    i = 0
    homopolymers_coords = []
    homopolymers_lengthes = []
    #print(seq_length)
    prev_end = 0
    while i < seq_length:
        if seq[i] == nucleotide:
            end, new_start = find_homopolymer_end(seq, nucleotide, seq_length, i, search_type=search_type,
                                                  max_single_insert_size=max_single_insert_size,
                                                  max_total_insert_length=max_total_insert_length,
                                                  max_number_of_insertions=max_number_of_insertions)
            # print(end, i)
            length = end - i
            if homopolymers_coords:
                prev_end = homopolymers_coords[-1][1]

            if length >= min_size and end > prev_end:
                homopolymers_coords.append((i, end))
                homopolymers_lengthes.append(length)

            i = end
            if new_start is not None:
                i = new_start
            continue
        i += 1

    homopolymers_lengthes = np.array(homopolymers_lengthes)
    return homopolymers_coords, homopolymers_lengthes

# ----------------------Filters-----------------------


def filter_sequences(record_dict, expression):
    record_id_list = []
    for record_id in record_dict:
        if expression(record_dict[record_id]):
            record_id_list.append(record_id)
    return record_id_list

# --------------------Generators----------------------


def rev_com_generator(record_dict, yield_original_record=False, rc_suffix="_rc", add_rc_suffix="False"):
    for record_id in record_dict:
        reverse_seq = record_dict[record_id].seq.reverse_complement()
        record = deepcopy(record_dict[record_id])
        record.seq = reverse_seq
        record.id = record_id + rc_suffix if yield_original_record or add_rc_suffix else record.id
        if yield_original_record:
            yield record_dict[record_id]
        yield record


def find_gaps(record_dict):
    gap_reg_exp = re.compile("N+", re.IGNORECASE)
    gaps_dict = {}
    for region in record_dict:
        gaps_dict[region] = SeqRecord(seq=record_dict[region].seq,
                                      id=record_dict[region].id,
                                      description=record_dict[region].description)
        gaps = gap_reg_exp.finditer(str(record_dict[region].seq))  # iterator with
        for match in gaps:
            gaps_dict[region].features.append(SeqFeature(FeatureLocation(match.start(), match.end()),
                                              type="gap", strand=None))
    return gaps_dict

# legacy function, moved to class SequenceRoutines as method, possibly all usages were changed
"""

def record_by_id_generator(record_dict, id_list):
    for record_id in id_list:
        if record_id in record_dict:
            yield record_dict[record_id]
        else:
            sys.stderr.write("Not found: %s\n" % record_id)
"""


def record_by_expression_generator(record_dict, expression=None, id_file="passed_records.ids"):
    """
    :param record_dict: dictionary containing Biopython SeqRecords as values
    :param expression: function to apply to all records in record_dict. If it returns True record will be yielded
    :param id_file: file to write ids of records passed expression
    :return: None
    """
    id_fd = open(id_file, "w")
    if expression is None:
        for record_id in record_dict:
            id_fd.write(record_id + "\n")
            yield record_dict[record_id]
    else:
        for record_id in record_dict:
            if expression(record_dict[record_id]):
                id_fd.write(record_id + "\n")
                yield record_dict[record_id]
    id_fd.close()


def record_generator(annotations_dict, sequence_dict, feature_types_list):
    for record_id in annotations_dict:
        for feature in annotations_dict[record_id].features:
            if feature.type in feature_types_list:
                sequence = feature.extract(sequence_dict[record_id].seq)
                #record = SeqRecord(sequence, id=feature.id)
                #print(record)
                yield SeqRecord(sequence, id=feature.id, description=feature.qualifiers["Name"][0] \
                      if "Name" in feature.qualifiers else "")


def get_kmer_dict(sequence, kmer_length, start=1, end=None):  # positions are one-based
    kmer_dict = OrderedDict()
    length = len(sequence)
    if end > length:
        raise ValueError("End position should be less than length")
    stop = length if end is None else end
    for i in range(start-1, stop-kmer_length):
        kmer_dict[i+1] = sequence[i:i+kmer_length]
    return kmer_dict


def get_kmer_dict_as_seq_records(sequence, kmer_length, start=1, end=None, id_prefix="id"):  # positions are one-based
    kmer_dict = OrderedDict()
    length = len(sequence)
    if end > length:
        raise ValueError("End position should be less than length")
    stop = length if end is None else end
    for i in range(start-1, stop-kmer_length):
        record_id = "%s_%i-%i" % (id_prefix, i + 1, i + kmer_length)
        kmer_dict[record_id] = SeqRecord(seq=sequence[i:i+kmer_length], id=record_id)
    return kmer_dict


def split_record_ids_by_expression(sequence_dict, expression):
    """
    splits record ids based on parameter of record
    """
    id_dict = OrderedDict()
    for record_id in sequence_dict:

        value = expression(sequence_dict[record_id])
        if value not in id_dict:
            id_dict[value] = [record_id]
        else:
            id_dict[value].append(record_id)

    return id_dict

if __name__ == "__main__":
    sequence = "CTGGCAAAGACCCAAACATCGACCACATCGAACAGCCACACCACCACCAACACTGTGCACCACTTCCGATTTCCAGCACCCCTTTTTGCCACTCTTTTTACGTAGTTTTGGCCATGCCTAGTTGTTTCCCAGTAGTCAACTTAAACGTATTTATTTTAATAAATTTCCACAAGGTTC"
    coords, length = find_homopolymers(sequence, "T", min_size=2, search_type="non_perfect",
                      max_single_insert_size=1, max_total_insert_length=None, max_number_of_insertions=2)
    print(coords)
    print(length)
    for coord in coords:
        print(sequence[coord[0]:coord[1]])