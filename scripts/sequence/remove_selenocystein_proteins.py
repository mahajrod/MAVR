#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
import os

from Bio import SeqIO

from Routines.Sequence import record_by_id_generator
from Routines.Sequence import record_by_expression_generator
parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input_file", required=True,
                    help="Input file with sequences")
parser.add_argument("-o", "--output_file", action="store", dest="output_file", required=True,
                    help="Output file with proteins without selenocystein")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of input and output files. Allowed formats genbank, fasta(default)")

args = parser.parse_args()

tmp_index_file = "temp.idx"

print("Parsing %s..." % args.input_file)
sequence_dict = SeqIO.index_db(tmp_index_file, args.input_file, format=args.format)
"""
all_ids = set(sequence_dict.keys())
selenocystein_ids = set()
for record_id in sequence_dict:
    if "U" in sequence_dict[record_id].seq:
        selenocystein_ids.add(record_id)
no_selenocystein_ids = all_ids - selenocystein_ids

record_by_expression_generator(record_dict, expression)
"""
SeqIO.write(record_by_expression_generator(sequence_dict, lambda record: "U" not in record.seq),
            args.output_file, args.format)
os.remove(tmp_index_file)



