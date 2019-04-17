#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import argparse
from multiprocessing import Pool
from Bio import SearchIO
from RouToolPa.Tools.LinuxTools import CGAS
from RouToolPa.Routines import FileRoutines # make_list_of_path_to_files, split_filename, check_path, save_mkdir



def make_list_of_path_to_files_from_comma_sep_string(string):
    return FileRoutines.make_list_of_path_to_files(string.split(","))

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    type=make_list_of_path_to_files_from_comma_sep_string,
                    help="Comma-separated list of files with hmm output")
parser.add_argument("-f", "--format", action="store", dest="format", required=True,
                    help="Format of input hmm file.")
parser.add_argument("-o", "--output_file", action="store", dest="output", default="stdout",
                    help="Output file")
parser.add_argument("-t", "--threads", action="store", dest="threads", default=1, type=int,
                    help="Number of threads to handle input")
parser.add_argument("-d", "--top_hits_dir", action="store", dest="top_hits_dir", default="top_hits_dir/",
                    type=FileRoutines.check_path,
                    help="Directory to write intermediate(splited) output")
parser.add_argument("-r", "--retain_splited_output", action="store_true", dest="retain",
                    help="Retain splited output")

args = parser.parse_args()

FileRoutines.safe_mkdir(args.top_hits_dir)


def handle_input(filename):
    sys.stdout.write("Handling %s\n" % filename)
    prefix = FileRoutines.split_filename(filename)[1]
    index_file = "%s.tmp.idx" % prefix
    hmm_dict = SearchIO.index_db(index_file, filename, args.format)
    if args.output == "stdout":
        out_fd = sys.stdout
    else:
        out_fd = open("%s%s.top_hits" % (args.top_hits_dir, prefix), "w")
        out_fd.write("#query\thit\tevalue\tbitscore\n")
    for family in hmm_dict:
        #print hmm_dict[key]
        for hit in hmm_dict[family]:
            if hit.is_included:
                out_fd.write("%s\t%s\t%s\t%s\n" % (family, hit.id, hit.evalue, hit.bitscore))
    if args.output != "stdout":
        out_fd.close()

    os.remove(index_file)

if args.output == "stdout":
    sys.stdout.write("#query\thit\tevalue\tbitscore\n")

process_pool = Pool(args.threads)
process_pool.map(handle_input, args.input)

if args.output != "stdout":
    CGAS.cat(["%s%s" % (args.top_hits_dir, filename) for filename in os.listdir(args.top_hits_dir)], output=args.output)

if not args.retain:
    os.remove(args.top_hits_dir)

"""
hmm_dict = SearchIO.index_db("temp.idx", args.input, args.format)
out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")
out_fd.write("#query\thit\tevalue\tbitscore\n")
for family in hmm_dict:
    #print hmm_dict[key]
    for hit in hmm_dict[family]:
        if hit.is_included:
            out_fd.write("%s\t%s\t%s\t%s\n" % (family, hit.id, hit.evalue, hit.bitscore))
if args.output != "stdout":
    out_fd.close()

os.remove("temp.idx")
"""