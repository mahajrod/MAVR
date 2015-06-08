#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
import os

from Bio import SeqIO

from Routines.Sequence import record_by_id_generator
parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input_file",
                    help="Input file with sequences")
parser.add_argument("-o", "--output_file_prefix", action="store", dest="out_prefix",
                    default="selenocystein_proteins",
                    help="Prefix_of_output_files")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of input and output files. Allowed formats genbank, fasta(default)")


args = parser.parse_args()

tmp_index_file = "temp.idx"

print("Parsing %s..." % args.input_file)
sequence_dict = SeqIO.index_db(tmp_index_file, args.input_file, format=args.format)
selenocystein_ids = []
with open(args.out_prefix + "_id.t", "w") as out_fd:
    for record_id in sequence_dict:
        if "U" in sequence_dict[record_id].seq:
            selenocystein_ids.append(record_id)
            out_fd.write(record_id + "\n")
SeqIO.write(record_by_id_generator(sequence_dict, selenocystein_ids), "%s_seq.fasta" % args.out_prefix, args.format)
os.remove(tmp_index_file)



