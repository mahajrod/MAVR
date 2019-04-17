#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import argparse
from Bio import SeqIO
from RouToolPa.Routines.File import make_list_of_path_to_files
from RouToolPa.Routines.Sequence import record_by_expression_generator



def check_if_protein_is_complete(record):

    splited_description = record.description.split()
    complete_orf = None
    for entry in splited_description:
        if ("type" in entry) and (":" in entry):
            complete_orf = entry.split(":")[1] == "complete"

    return complete_orf

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file_list", action="store", dest="input", required=True,
                    type=lambda s: make_list_of_path_to_files(s.split(",")),
                    help="Comma-separated list of input fasta files or directories with"
                         "peptide output from TransDecoder")
parser.add_argument("-o", "--output_file", action="store", dest="output", default="stdout",
                    help="Output fasta file with complete peptides")

parser.add_argument("-d", "--id_file", action="store", dest="id_file", default="complete.pep.ids",
                    help="File to write ids of complete peptides. Default - complete.pep.ids")

args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

tmp_index_file = "temp.idx"

print("Parsing %s..." % (args.input if isinstance(args.input, str) else ",".join(args.input)))
sequence_dict = SeqIO.index_db(tmp_index_file, args.input, format="fasta")

SeqIO.write(record_by_expression_generator(sequence_dict, expression=check_if_protein_is_complete,
                                           id_file=args.id_file),
            out_fd, format="fasta")

os.remove(tmp_index_file)

if args.output != "stdout":
    out_fd.close()


