#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import argparse
from multiprocessing import Pool

from CustomCollections.GeneralCollections import IdList
from Routines.File import check_path

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", default="stdin",
                    help="Input file with ids of families. Default: stdin")
parser.add_argument("-t", "--threads", action="store", dest="threads", default=10, type=int,
                    help="Number of threads to download. Default: 10")
parser.add_argument("-o", "--output_directory", action="store", dest="output_dir", default="fam_data/",
                    help="Directory to output files.")

parser.add_argument("-l", "--alignment", action="store_true", dest="alignment", default=False,
                    help="Download alignments")
parser.add_argument("-r", "--tree", action="store_true", dest="tree", default=False,
                    help="Download trees")
parser.add_argument("-h", "--hmm", action="store_true", dest="hmm", default=False,
                    help="Download hmms")
parser.add_argument("-a", "--alignment", action="store_true", dest="all", default=False,
                    help="Download all: alignment, tree, hmm")
args = parser.parse_args()

args.output_dir = check_path(args.output_dir)
if (not args.alignment) and (not args.tree) and (not args.hmm):
    args.all = True

in_fd = sys.stdin if args.input == "stdin" else open(args.input, "r")

family_ids = IdList()
family_ids.read(in_fd)
if args.input != "stdin":
    in_fd.close()


def download_data(fam_id):
    print("Downloading %s family" % fam_id)
    alignment_options = "http://www.treefam.org/family/%s/alignment -o %s%s.fasta" % \
                        (fam_id, args.output_dir, fam_id)
    tree_options = "http://www.treefam.org/family/%s/tree/newick -o %s%s.nwk" % \
                   (fam_id, args.output_dir, fam_id)
    hmm_options = "http://www.treefam.org/family/%s/hmm -o %s%s.hmm" % \
                  (fam_id, args.output_dir, fam_id)

    if args.all or args.alignment:
        os.system("wget %s" % alignment_options)
    if args.all or args.tree:
        os.system("wget %s" % tree_options)
    if args.all or args.hmm:
        os.system("wget %s" % hmm_options)


pool = Pool(args.threads)
pool.map(download_data, family_ids)

