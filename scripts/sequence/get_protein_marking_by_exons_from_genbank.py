#!/usr/bin/env python
__author__ = 'mahajrod'
import os
import sys
import argparse

from Bio import SeqIO

from Routines import SequenceRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", type=SequenceRoutines.make_list_of_path_to_files,
                    help="Comma-separated list of genbank files/directories with annotations")
parser.add_argument("--fast_parsing", action="store_true", dest="fast_parsing",
                    help="Fast parsing mode - high memory consumption. Default: false")
parser.add_argument("-o", "--output", action="store", dest="output",
                    default="out.t", help="output file - default: out.t.")

args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

record_dict = SeqIO.to_dict(SeqIO.parse(args.input, format="genbank")) if args.fast_parsing else SeqIO.index_db("temp_index.idx", args.input, format="genbank")

SequenceRoutines.get_protein_marking_by_exons_from_genbank(record_dict, args.output,
                                                           protein_id_field_in_cds_feature="protein_id")

#os.remove("temp_index.idx")
