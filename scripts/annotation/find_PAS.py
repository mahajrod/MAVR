#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse
from Bio import SeqIO
from RouToolPa.Routines.File import make_list_of_path_to_files


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    type=lambda s: make_list_of_path_to_files(s.split(",")),
                    help="Comma-separated list of fasta files with sequences")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output bed file with coordinates")

args = parser.parse_args()

PAS_sites = [
             "AATAAA",
             "ATTAAA",
             "TATAAA",
             "AGTAAA",
             "AAGAAA",
             "GATAAA",
             "CATAAA",
             "AATATA",
             "AATACA",
             "AATGAA",
             "ACTAAA",
             "AATAGA"]

sequence_dict = SeqIO.index_db("temp.idx", args.input, format="fasta")

with open(args.output, "w") as out_fd:
    for record_id in sequence_dict:
        record_length = len(sequence_dict[record_id])
        seq = sequence_dict[record_id].seq.upper()
        for i in range(0, record_length - 6):
            if seq[i:i+6] in PAS_sites:
                out_fd.write("%s\t%i\t%i\t%s\n" % (record_id, i+1, i+6, seq[i:i+6]))

os.remove("temp.idx")
