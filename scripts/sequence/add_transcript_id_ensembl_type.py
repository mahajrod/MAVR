#!/usr/bin/env python
__author__ = 'mahajrod'
import os
import sys
import argparse

from RouToolPa.Collections.General import SynDict
from RouToolPa.Parsers.Sequence import CollectionSequence


parser = argparse.ArgumentParser()

parser.add_argument("-f", "--format", action="store", dest="format",
                    default="fasta", help="format of file with sequences - default: fasta.")
parser.add_argument("-i", "--input", action="store", dest="input",
                    help="File with sequences")
parser.add_argument("-o", "--output", action="store", dest="output",
                    default="stdout", help="output file")
parser.add_argument("-s", "--scaffold_syn_file", action="store", dest="scaffold_syn_file",
                    help="File with scaffold id synonyms")
parser.add_argument("-k", "--syn_file_key_column", action="store", dest="syn_file_key_column",
                    default=0, type=int,
                    help="Column(0-based) with key(current id) for scaffolds in synonym file. Default: 0")
parser.add_argument("-v", "--syn_file_value_column", action="store", dest="syn_file_value_column",
                    default=1, type=int,
                    help="Column(0-based) with value(synonym id) for scaffolds in synonym file synonym. Default: 1")

args = parser.parse_args()

fasta_collection = CollectionSequence(in_file=args.input, format=args.format)

syn_dict = SynDict(filename=args.scaffold_syn_file,
                   key_index=args.syn_file_key_column,
                   value_index=args.syn_file_value_column)

for seq_id in fasta_collection.records:
    if seq_id in syn_dict:
        fasta_collection.description[seq_id] += "{0}transcript:{1}".format(" " if fasta_collection.description[seq_id] != "" else "",
                                                                           syn_dict[seq_id])
fasta_collection.write(args.output)
