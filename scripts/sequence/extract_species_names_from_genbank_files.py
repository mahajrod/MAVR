#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import argparse
from Bio import SeqIO
from RouToolPa.Routines.File import make_list_of_path_to_files



parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file_list", action="store", dest="input", required=True,
                    type=lambda s: make_list_of_path_to_files(s.split(",")),
                    help="Comma-separated list of genbank files/directories with sequences")
parser.add_argument("-o", "--output_file", action="store", dest="output", default="stdout",
                    help="Output file with species")

args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

tmp_index_file = "temp.idx"

sequence_dict = SeqIO.index_db(tmp_index_file, args.input, format="genbank")

for record_id in sequence_dict:
    #print sequence_dict[record_id].annotations
    out_fd.write("%s\t%s\t\n" % (record_id, sequence_dict[record_id].annotations['organism']))

os.remove(tmp_index_file)

if args.output != "stdout":
    out_fd.close()

