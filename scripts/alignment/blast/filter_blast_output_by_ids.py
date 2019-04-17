#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from Bio import SearchIO
from Bio.SearchIO import QueryResult
from RouToolPa.Routines.File import read_ids


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", default="stdin",
                    help="Input file with BLAST results. Default: stdin")
parser.add_argument("-o", "--output_file", action="store", dest="output", default="stdout",
                    help="Output file with filtered BLAST results. Default: stdout")
parser.add_argument("-f", "--format", action="store", dest="format", default="blast-tab",
                    help="Format of BLAST file. Allowed: blast-tab, blast-xml, blat-psl, hmmer3-tab, hmmer3-domtab. Default: blast-tab")
parser.add_argument("-m", "--mode", action="store", dest="mode", default="both",
                    help="Operating mode of filtering: query - filter query by ids, database - filter database sequences by ids, both - filter both. Default: both ")
parser.add_argument("-w", "--white_list_id_file", action="store", dest="white_list_id_file",
                    help="File with ids from white list. Cannot be set with -b/--black_list_id_file")
parser.add_argument("-b", "--black_list_id_file", action="store", dest="black_list_id_file",
                    help="File with ids from black list. Cannot be set with -w/--white_list_id_file.")
args = parser.parse_args()

in_fd = sys.stdin if args.input == "stdin" else open(args.input, "r")
out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

allowed_formats = ["blast-tab",
                   "blast-xml",
                   "blat-psl",
                   "hmmer3-tab",
                   "hmmer3-domtab"]
allowed_modes = ["query",
                 "database",
                 "both"]

if args.format not in allowed_formats:
    raise ValueError("This format of input file is not allowed")
elif (args.white_list_id_file is None) and (args.black_list_id_file is None):
    raise ValueError("Both files with ids from black and white list were not set")
elif (args.white_list_id_file is not None) and (args.black_list_id_file is not None):
    raise ValueError("Both files with ids from black and white list were set")
elif args.mode not in allowed_modes:
    raise ValueError("This filtering mode is not allowed")

if args.white_list_id_file:
    white_list = read_ids(args.white_list_id_file)
    if args.mode == "query":
        def iterator(blast_dict):
            for entry in blast_dict:
                if entry in white_list:
                    yield blast_dict[entry]
    elif args.mode == "database":
        def iterator(blast_dict):
            for entry in blast_dict:
                entry_hits = []
                for hit in blast_dict[entry].hits:
                    if hit.id in white_list:
                        # filter hits
                        entry_hits.append(hit)
                if entry_hits:
                    yield QueryResult(hits=entry_hits, id=entry)

    elif args.mode == "both":
        def iterator(blast_dict):
            for entry in blast_dict:
                if entry in white_list:
                    entry_hits = []
                    for hit in blast_dict[entry].hits:
                        if hit.id in white_list:
                            # filter hits
                            entry_hits.append(hit)
                    if entry_hits:
                        yield QueryResult(hits=entry_hits, id=entry)
else:
    black_list = read_ids(args.black_list_id_file)
    if args.mode == "query":
        def iterator(blast_dict):
            for entry in blast_dict:
                if entry not in black_list:
                    yield blast_dict[entry]
    elif args.mode == "database":
        def iterator(blast_dict):
            for entry in blast_dict:
                entry_hits = []
                for hit in blast_dict[entry].hits:
                    if hit.id not in black_list:
                        # filter hits
                        entry_hits.append(hit)
                if entry_hits:
                    yield QueryResult(hits=entry_hits, id=entry)

    elif args.mode == "both":
        def iterator(blast_dict):
            for entry in blast_dict:
                if entry not in black_list:
                    entry_hits = []
                    for hit in blast_dict[entry].hits:
                        if hit.id not in black_list:
                            # filter hits
                            entry_hits.append(hit)
                    if entry_hits:
                        yield QueryResult(hits=entry_hits, id=entry)

blast_results = SearchIO.index(args.input, args.format)

SearchIO.write(iterator(blast_results), args.output, args.format)
if args.output != "output":
    out_fd.close()
if args.input != "stdin":
    in_fd.close()