#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os
import argparse
from collections import OrderedDict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from RouToolPa.Routines.File import make_list_of_path_to_files
from RouToolPa.Routines.Sequence import record_by_expression_generator


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    type=lambda s: make_list_of_path_to_files(s.split(",")),
                    help="Comma-separated list files/directories with sequences")

parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-e", "--allow_overlaps", action="store_true", dest="allow_overlaps",
                    help="Allow overlaps of terminal fragments")
parser.add_argument("-l", "--left_seq_size", action="store", dest="left_seq_size", required=True, type=int,
                    help="Size of left terminal fragment")
parser.add_argument("-r", "--right_seq_size", action="store", dest="right_seq_size", required=True, type=int,
                    help="Size of right terminal fragment")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of input files. Allowed: fasta, genbank. Default: fasta")

args = parser.parse_args()

print("Parsing %s" % args.input)
sequence_dict = SeqIO.index_db("temp.idx", args.input, format=args.format)

left_fragments_dict = OrderedDict()
right_fragments_dict = OrderedDict()
overlap_counter = 0
number_of_records = len(sequence_dict)
for record_id in sequence_dict:
    record_size = len(sequence_dict[record_id].seq)

    if (record_size < (args.left_seq_size + args.right_seq_size)) and (not args.allow_overlaps):
        print("%s - terminal fragments overlap" % record_id)
        overlap_counter += 1
        continue

    #left_id = "%s_1-%i" % (record_id, args.left_seq_size)
    #right_id = "%s_%i-%i" % (record_id, record_size-args.right_seq_size+1, record_size)

    left_record = SeqRecord(seq=sequence_dict[record_id].seq[:args.left_seq_size], id=record_id) # id=left_id)
    right_record = SeqRecord(seq=sequence_dict[record_id].seq[-args.right_seq_size:], id=record_id) # id=right_id)

    left_fragments_dict[record_id] = left_record
    right_fragments_dict[record_id] = right_record

SeqIO.write(record_by_expression_generator(left_fragments_dict, lambda x: True),
            "%s_left_fragments.fasta" % args.output_prefix, format="fasta")
SeqIO.write(record_by_expression_generator(right_fragments_dict, lambda x: True),
            "%s_right_fragments.fasta" % args.output_prefix, format="fasta")

print("Totally %i/%i records with overlaps" % (overlap_counter, number_of_records))
os.remove("temp.idx")
