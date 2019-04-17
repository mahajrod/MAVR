#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse
from Bio import SeqIO
from RouToolPa.Routines import MtDNARoutines, FileRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    type=lambda s: FileRoutines.make_list_of_path_to_files(s.split(",")),
                    help="Comma separated list of genbank files/directories")
#parser.add_argument("-o", "--output_file_prefix", action="store", dest="output_prefix", required=True,
#                    help="Prefix of output files")
parser.add_argument("-f", "--format", action="store", dest="format", default="genbank",
                    help="Format of input and output file. Allowed formats genbank(default), fasta")

args = parser.parse_args()

record_dict = SeqIO.index_db("tmp.idx", args.input, format=args.format)

MtDNARoutines.split_mitochondrion_genome_by_genes(record_dict, black_list=[])

os.remove("tmp.idx")

