__author__ = 'mahajrod'

from copy import deepcopy
from collections import OrderedDict

from Routines.Functions import output_dict


def get_lengths(record_dict, out_file="lengthes.t", write=False):
    lengths_dict = OrderedDict({})
    for record_id in record_dict:
        lengths_dict[record_id] = len(record_dict[record_id])

    if write:
        output_dict(lengths_dict, out_file=out_file, write=write)

    return lengths_dict

#--------------------Generators----------------------


def rev_com_generator(record_dict):
    for record_id in record_dict:
        reverse_seq = record_dict[record_id].seq.reverse_complement()
        record = deepcopy(record_dict[record_id])
        record.seq = reverse_seq
        yield record