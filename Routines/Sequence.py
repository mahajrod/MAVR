__author__ = 'mahajrod'

from copy import deepcopy
from collections import OrderedDict

import numpy as np

from Bio.SeqRecord import SeqRecord

from CustomCollections.GeneralCollections import TwoLvlDict
from Routines.Functions import output_dict


def get_lengths(record_dict, out_file="lengths.t", write=False):
    lengths_dict = OrderedDict({})
    for record_id in record_dict:
        lengths_dict[record_id] = len(record_dict[record_id])

    if write:
        output_dict(lengths_dict, out_file=out_file, write=write)

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


#--------------------Search--------------------------
#              !!!!IMPORTANT!!!!!!
#in this section python notation for coordinates inside sequence is used
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


def record_by_id_generator(record_dict, id_list):
    for record_id in id_list:
        if record_id in record_dict:
            yield record_dict[record_id]
        else:
            print (record_id)


def record_generator(annotations_dict, sequence_dict, feature_types_list):
    for record_id in annotations_dict:
        for feature in annotations_dict[record_id].features:
            if feature.type in feature_types_list:
                sequence = feature.extract(sequence_dict[record_id].seq)
                #record = SeqRecord(sequence, id=feature.id)
                #print(record)
                yield SeqRecord(sequence, id=feature.id, description=feature.qualifiers["Name"][0] \
                      if "Name" in feature.qualifiers else "")

if __name__ == "__main__":
    sequence = "CTGGCAAAGACCCAAACATCGACCACATCGAACAGCCACACCACCACCAACACTGTGCACCACTTCCGATTTCCAGCACCCCTTTTTGCCACTCTTTTTACGTAGTTTTGGCCATGCCTAGTTGTTTCCCAGTAGTCAACTTAAACGTATTTATTTTAATAAATTTCCACAAGGTTC"
    coords, length = find_homopolymers(sequence, "T", min_size=2, search_type="non_perfect",
                      max_single_insert_size=1, max_total_insert_length=None, max_number_of_insertions=2)
    print(coords)
    print(length)
    for coord in coords:
        print(sequence[coord[0]:coord[1]])