#!/usr/bin/env python
__author__ = 'mahajrod'
import os
import argparse

from SeqAnalysis.SeqAnalysis import *

parser = argparse.ArgumentParser()
parser.add_argument('-i', action='store', dest='input_file', help='Input file with sequences in .gb format')

args = parser.parse_args()
workdir = os.getcwd()

if args.input_file[-3:] == ".gb":
	input_prefix = args.input_file[:-3]
elif args.input_file[-8:] == ".genbank":
	input_prefix = args.input_file[:-8]
else:
	raise ValueError("Input filetype was not recognized. Only .gb or .genbank extentions are allowed")

input_file = args.input_file
input_index = input_prefix + ".idx"
filtered_prefix = input_prefix.split("/")[-1] + "_filtered"

filter_taxa_list = ["unclassified", " x ", "sp.", "artificial sequences", "synthetic construct"]

print("Parsing %s..." % input_file)

record_dict, taxonomy_dict = get_statistics([input_file],
                                            input_index,
                                            taxonomy_file_prefix=input_prefix)



#count_species(record_dict, output_filename=input_prefix + "_species.count")
filter_by_taxa(record_dict,
               filter_taxa_list,
               filter_type="black_list",
               output_filename=filtered_prefix+".gb",
               output_type="genbank")

#filtered_record_dict = SeqIO.index_db(filtered_prefix + ".idx", [filtered_prefix + ".gb"], "genbank")

filtered_record_dict, taxonomy_dict = get_statistics([filtered_prefix + ".gb"],
                                            filtered_prefix + ".idx",
                                            taxonomy_file_prefix=filtered_prefix)

count_species(filtered_record_dict, output_filename=filtered_prefix+"_species.count")
get_taxonomy(filtered_record_dict, output_prefix=filtered_prefix)

sort_by_number(filtered_prefix+"_species.count",
               filtered_prefix,
               counts_list=[1, 5, 10, 15, 25])


"""
filter_by_taxa(filtered_record_dict,
               []
               filter_type="white_list",
               output_filename="Pisces_filtered_by_taxa.gb",
               filtered_out_filename="non_pisces_by_taxa.gb",
               output_type="genbank",
               return_record_generator=False)
"""