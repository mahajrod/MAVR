#!/usr/bin/env python
__author__ = 'mahajrod'
import os
import argparse

from Bio import SeqIO

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input",
                    help="fasta file with sequences")
parser.add_argument("-e", "--extension", action="store", dest="ext", default=".fasta",
                    help="Extension of output fasta files. Default: .fasta")
args = parser.parse_args()

record_dict = SeqIO.index_db("temp_index.idx", [args.input], format="fasta")

for record_id in record_dict:
    SeqIO.write([record_dict[record_id]], record_id + args.ext, "fasta")

os.remove("temp_index.idx")
