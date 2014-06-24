#!/usr/bin/env python

import os
import argparse

from SeqAnalysis.SeqAnalysis import *

parser = argparse.ArgumentParser()
parser.add_argument('-t', action='store', dest='taxa_list', help='Group of taxa to be retained and splited.Taxa inside group should be comma-separated')
parser.add_argument('-i', action='store', dest='input_file', help='Input file with sequences in .gb format')

args = parser.parse_args()

taxa_list = args.taxa_list.split(",")

if args.input_file[-3:] == ".gb":
	input_prefix = args.input_file[:-3]
elif args.input_file[-8:] == ".genbank":
	input_prefix = args.input_file[:-8]
else:
	raise ValueError("Input filetype was not recognized. Only .gb or .genbank extentions are allowed")

input_file = args.input_file
input_index = input_prefix + ".idx"

output_suffix = "_complete_genome"
print(input_index )
print(input_file)
record_dict = SeqIO.index_db(input_index, [input_file], "genbank")
os.system("mkdir -p splited")
os.chdir("splited")

split_by_taxa(record_dict, taxa_list, output_suffix, output_type="genbank")
	
