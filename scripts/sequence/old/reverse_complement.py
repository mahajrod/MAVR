#!/usr/bin/env python
__author__ = 'mahajrod'
import os
import argparse
from Bio import SeqIO
from RouToolPa.Routines.Sequence import rev_com_generator



parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input",
                    help="fasta file with sequences")
parser.add_argument("-w", "--write_original", action="store_true", dest="write_original", default=False,
                    help="Write original records")
parser.add_argument("-o", "--output", action="store", dest="output",
                    help="fasta file with reverse complement sequences")

args = parser.parse_args()

record_dict = SeqIO.index_db("temp_index.idx", [args.input], format="fasta")

SeqIO.write(rev_com_generator(record_dict, yield_original_record=args.write_original), args.output, "fasta")

os.remove("temp_index.idx")
