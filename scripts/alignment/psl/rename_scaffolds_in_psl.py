#!/usr/bin/env python
__author__ = 'mahajrod'
import sys
import argparse

import pandas as pd
from RouToolPa.GeneralRoutines import FileRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-p", "--psl", action="store", dest="psl", default=sys.stdin,
                    help="Input psl file. Default: stdin")

parser.add_argument("-s", "--query_syn_file", action="store", dest="query_syn_file",
                    help="File with query scaffold id synonyms")
parser.add_argument("--query_syn_file_key_column", action="store", dest="query_syn_file_key_column",
                    default=0, type=int,
                    help="Column(0-based) with key(current id) for scaffolds in query synonym file")
parser.add_argument("--query_syn_file_value_column", action="store", dest="query_syn_file_value_column",
                    default=1, type=int,
                    help="Column(0-based) with value(synonym id) for scaffolds in query synonym file synonym")

parser.add_argument("-y", "--target_syn_file", action="store", dest="target_syn_file",
                    help="File with target scaffold id synonyms")
parser.add_argument("--target_syn_file_key_column", action="store", dest="target_syn_file_key_column",
                    default=0, type=int,
                    help="Column(0-based) with key(current id) for scaffolds in target synonym file")
parser.add_argument("--target_syn_file_value_column", action="store", dest="target_syn_file_value_column",
                    default=1, type=int,
                    help="Column(0-based) with value(synonym id) for scaffolds in target synonym file synonym")

parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output file")

args = parser.parse_args()

QUERY_SCAFFOLD_ID_COLUMN = 9
TARGET_SCAFFOLD_ID_COLUMN = 13

if args.query_syn_file:
    query_syn_df = pd.read_csv(args.query_syn_file, sep="\t", usecols=(args.query_syn_file_key_column,
                                                                       args.query_syn_file_value_column),
                               index_col=args.query_syn_file_key_column,  header=None)
    query_syn_df.columns = ["syn"]
    query_syn_df.index.name = "key"

if args.target_syn_file:
    target_syn_df = pd.read_csv(args.target_syn_file, sep="\t", usecols=(args.target_syn_file_key_column,
                                                                         args.target_syn_file_value_column),
                                index_col=args.target_syn_file_key_column, header=None)
    target_syn_df.columns = ["syn"]
    target_syn_df.index.name = "key"

#print(query_syn_df)

if args.query_syn_file and args.target_syn_file:
    def rename_fuction(line_list):
        if line_list[9] in query_syn_df.index:
            line_list[9] = query_syn_df.loc[line_list[9], "syn"]
        if line_list[13] in target_syn_df.index:
            line_list[13] = target_syn_df.loc[line_list[13], "syn"]
elif args.query_syn_file:
    def rename_fuction(line_list):
        if line_list[9] in query_syn_df.index:
            line_list[9] = query_syn_df.loc[line_list[9], "syn"]
elif args.target_syn_file:
    def rename_fuction(line_list):
        if line_list[13] in target_syn_df.index:
            line_list[13] = target_syn_df.loc[line_list[13], "syn"]
else:
    def rename_fuction(line_list):
        return None

with FileRoutines.metaopen(args.psl, "r") as in_fd, FileRoutines.metaopen(args.output, "w") as out_fd:
    for line in in_fd:
        if line[0] == "#":
            out_fd.write(line)
            continue
        line_list = line.split("\t")
        rename_fuction(line_list)
        out_fd.write("\t".join(line_list))
