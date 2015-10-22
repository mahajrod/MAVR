#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import argparse

from collections import OrderedDict

from Bio import SeqIO
from BCBio import GFF

from Routines.File import read_ids, make_list_of_path_to_files
from Routines.Sequence import record_by_id_generator


def iterate_gff_parsing_generator(gff_parsing_generator, id_list):
    for record in gff_parsing_generator:
        if record.id in id_list:
            yield record

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file_list", action="store", dest="input", required=True,
                    type=lambda s: make_list_of_path_to_files(s.split(",")),
                    help="Comma-separated list of input files/directories with sequences")
parser.add_argument("-o", "--output_file", action="store", dest="output", default="stdout",
                    help="Output file with sequences")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of input and output files. Allowed formats genbank, fasta(default)")
parser.add_argument("-d", "--id_file", action="store", dest="id_file",
                    help="File with ids of sequences to extract")

args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

tmp_index_file = "temp.idx"

id_list = read_ids(args.id_file)

print("Parsing %s..." % (args.input if isinstance(args.input, str) else ",".join(args.input)))

sequence_dict = SeqIO.index_db(tmp_index_file, args.input, format=args.format)
SeqIO.write(record_by_id_generator(sequence_dict, id_list), out_fd, format=args.format)
os.remove(tmp_index_file)

if args.output != "stdout":
    out_fd.close()


