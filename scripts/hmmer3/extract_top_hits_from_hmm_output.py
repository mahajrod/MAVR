#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import argparse

from multiprocessing import Pool

from Bio import SearchIO

from CustomCollections.GeneralCollections import IdList
from Routines.File import make_list_of_path_to_files, split_filename, check_path, save_mkdir


def make_list_of_path_to_files_from_comma_sep_string(string):
    return make_list_of_path_to_files(string.split(","))

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_files", action="store", dest="input", required=True,
                    type=make_list_of_path_to_files_from_comma_sep_string,
                    help="Comma-separated list of directories and/or files containing hmm output")
parser.add_argument("-f", "--format", action="store", dest="format", required=True,
                    help="Format of input hmm file.")

parser.add_argument("-o", "--output_file", action="store", dest="output", default="stdout",
                    help="Output file")
parser.add_argument("-n", "--not_found", action="store", dest="not_found", default="not_found.ids",
                    help="File to write ids of not found queries")
parser.add_argument("-g", "--not_significant", action="store", dest="not_significant", default="not_significant.hits",
                    help="File to write ids of queries with not significant hits")
parser.add_argument("-t", "--threads", action="store", dest="threads", default=1, type=int,
                    help="Number of threads to handle input")
parser.add_argument("-d", "--top_hits_dir", action="store", dest="top_hits_dir", default="top_hits_dir/",
                    type=check_path,
                    help="Directory to write intermediate(splited) output")
parser.add_argument("-r", "--retain_splited_output", action="store_true", dest="retain",
                    help="Retain splited output")

args = parser.parse_args()

save_mkdir(args.top_hits_dir)


def handle_input(filename):
    sys.stdout.write("Handling %s\n" % filename)
    not_significant_ids = IdList()
    not_found_ids = IdList()

    prefix = split_filename(filename)[1]
    index_file = "%s.tmp.idx" % prefix
    hmm_dict = SearchIO.index_db(index_file, filename, args.format)
    if args.output == "stdout":
        out_fd = sys.stdout
    else:
        out_fd = open("%s%s.top_hits" % (args.top_hits_dir, prefix), "w")
        out_fd.write("#query\thit\tevalue\tbitscore\n")

    for query in hmm_dict:
        if hmm_dict[query].hits:
            if hmm_dict[query][0].is_included:
                out_fd.write("%s\t%s\t%s\t%s\n" % (query, hmm_dict[query][0].id, hmm_dict[query][0].evalue,
                                               hmm_dict[query][0].bitscore))
            else:
                not_significant_ids.append(query)
        else:
            not_found_ids.append(query)

    if args.output != "stdout":
        out_fd.close()

    os.remove(index_file)
    return not_significant_ids, not_found_ids

if args.output == "stdout":
    sys.stdout.write("#query\thit\tevalue\tbitscore\n")

process_pool = Pool(args.threads)
results = process_pool.map(handle_input, args.input)

if args.output != "stdout":
    with open(args.output, "w") as out_fd:
        out_fd.write("#query\thit\tevalue\tbitscore\n")
        for file_name in ["%s%s" % (args.top_hits_dir, filename) for filename in os.listdir(args.top_hits_dir)]:
            with open(file_name, "r") as file_fd:
                file_fd.readline()
                for line in file_fd:
                    out_fd.write(line)


nf_fd = open(args.not_found, "w")
ns_fd = open(args.not_significant, "w")

for not_significant, not_found in results:
    not_significant.write(ns_fd)
    not_found.write(nf_fd)

if not args.retain:
    os.remove(args.top_hits_dir)

nf_fd.close()
ns_fd.close()

"""
hmm_dict = SearchIO.index_db("temp.idx", args.input, args.format)
out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")
out_fd.write("#query\thit\tevalue\tbitscore\n")


for query in hmm_dict:
    #print hmm_dict[key]
    if hmm_dict[query].hits:
        if hmm_dict[query][0].is_included:
            out_fd.write("%s\t%s\t%s\t%s\n" % (query, hmm_dict[query][0].id, hmm_dict[query][0].evalue,
                                               hmm_dict[query][0].bitscore))
        else:
            ns_fd.write("%s\t%s\t%s\t%s\n" % (query, hmm_dict[query][0].id, hmm_dict[query][0].evalue,
                                               hmm_dict[query][0].bitscore))
    else:
        nf_fd.write("%s\n" % query)

if args.output != "stdout":
    out_fd.close()

os.remove("temp.idx")


"""