#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse
from Bio import SeqIO
from RouToolPa.Collections.General import SynDict, IdSet
from RouToolPa.Routines import FileRoutines


#from RouToolPa.Routines.File import read_ids, make_list_of_path_to_files, check_path, save_mkdir
from RouToolPa.Routines import SequenceRoutines


#def make_list_of_path_to_files_from_comma_sep_string(string):
#    return make_list_of_path_to_files(string.split(","))

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file_list", action="store", dest="input", required=True,
                    type=FileRoutines.make_list_of_path_to_files_from_string,
                    help="Comma-separated list of input files/directories with sequences")
parser.add_argument("-o", "--output_directory", action="store", dest="output", type=FileRoutines.check_path,
                    help="Directory to output groups_of sequences")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of input and output files. Allowed formats genbank, fasta(default)")
parser.add_argument("-e", "--extension", action="store", dest="extension",
                    help="Extension of output files. Default: equal to -f")
parser.add_argument("-d", "--id_file", action="store", dest="id_file",
                    help="File with groups of sequences to extract(.fam file).")

args = parser.parse_args()

FileRoutines.safe_mkdir(args.output)
args.extension = args.extension if args.extension else args.format
tmp_index_file = "temp.idx"

#id_list = read_ids(args.id_file)
id_list = IdSet(filename=args.id_file)

sequence_groups_id = SynDict()
sequence_groups_id.read(args.id_file, split_values=True)
#print("Parsing %s..." % args.input_file)
sequence_dict = SeqIO.index_db(tmp_index_file, args.input, format=args.format)
for group in sequence_groups_id:
    SeqIO.write(SequenceRoutines.record_by_id_generator(sequence_dict, sequence_groups_id[group],
                                                        verbose=True),
                "%s%s.%s" % (args.output, group, args.extension), format=args.format)

os.remove(tmp_index_file)
