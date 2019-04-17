#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import argparse
from multiprocessing import Pool
from RouToolPa.Collections.General import IdList
from RouToolPa.Routines import FileRoutines



parser = argparse.ArgumentParser()

parser.add_argument("-i", "--id_file", action="store", dest="input", default="stdin",
                    help="Input file with ids of families. Default: stdin")
parser.add_argument("-t", "--threads", action="store", dest="threads", default=10, type=int,
                    help="Number of threads to download. Default: 10")
parser.add_argument("-o", "--output_directory", action="store", dest="output_dir", default="fam_data/", type=FileRoutines.check_path,
                    help="Directory to output files.")
parser.add_argument("-l", "--alignment", action="store_true", dest="alignment", default=False,
                    help="Download alignments")
parser.add_argument("-r", "--tree", action="store_true", dest="tree", default=False,
                    help="Download trees")
parser.add_argument("-m", "--hmm", action="store_true", dest="hmm", default=False,
                    help="Download hmms")
parser.add_argument("-a", "--all", action="store_true", dest="all", default=False,
                    help="Download all: alignment, tree, hmm")
parser.add_argument("-s", "--store_logs", action="store_true", dest="store_logs", default=False,
                    help="Store download logs in directory set by -g/--logs_dir option")
parser.add_argument("-g", "--logs_dir", action="store", dest="logs_dir", default="logs", type=FileRoutines.check_path,
                    help="Directory with logs")
args = parser.parse_args()

FileRoutines.safe_mkdir(args.output_dir)
FileRoutines.safe_mkdir(args.logs_dir)

if (not args.alignment) and (not args.tree) and (not args.hmm):
    args.all = True

in_fd = sys.stdin if args.input == "stdin" else open(args.input, "r")

family_ids = IdList()
family_ids.read(in_fd)

if args.input != "stdin":
    in_fd.close()

absent_alignment_list = IdList()
absent_tree_list = IdList()
absent_hmm_list = IdList()


def download_data(fam_id):
    print("Downloading %s family" % fam_id)
    ali_log_file = "/dev/null" if not args.store_logs else "%s%s_alignment.log" % (args.logs_dir, fam_id)
    tree_log_file = "/dev/null" if not args.store_logs else "%s%s_tree.log" % (args.logs_dir, fam_id)
    hmm_log_file = "/dev/null" if not args.store_logs else "%s%s_hmm.log" % (args.logs_dir, fam_id)

    alignment_options = " -O %s%s.fasta -c -t 2000 -o %s http://www.treefam.org/family/%s/alignment " % \
                        (args.output_dir, fam_id, ali_log_file, fam_id)
    tree_options = " -O %s%s.nwk -c -t 2000 -o %s http://www.treefam.org/family/%s/tree/newick" % \
                   (args.output_dir, fam_id, tree_log_file, fam_id)
    hmm_options = "  -O %s%s.hmm -c -t 2000 -o %s http://www.treefam.org/family/%s/hmm" % \
                  (args.output_dir, fam_id, hmm_log_file, fam_id)

    if args.all or args.alignment:
        os.system("wget %s" % alignment_options)
    if args.all or args.tree:
        os.system("wget %s" % tree_options)
    if args.all or args.hmm:
        os.system("wget %s" % hmm_options)


pool = Pool(args.threads)
pool.map(download_data, family_ids)
pool.close()
for fam_id in family_ids:
    if args.all or args.alignment:
        if os.path.getsize("%s%s.fasta" % (args.output_dir, fam_id)) == 0:
            absent_alignment_list.append(fam_id)
    if args.all or args.tree:
        if os.path.getsize("%s%s.nwk" % (args.output_dir, fam_id)) == 0:
            absent_tree_list.append(fam_id)
    if args.all or args.hmm:
        if os.path.getsize("%s%s.hmm" % (args.output_dir, fam_id)) == 0:
            absent_hmm_list.append(fam_id)

if absent_alignment_list:
    absent_alignment_list.write("absent_alignments.ids")
    print("%i alignments were not downloaded" % len(absent_alignment_list))
if absent_tree_list:
    absent_tree_list.write("absent_trees.ids")
    print("%i trees were not downloaded" % len(absent_tree_list))
if absent_hmm_list:
    absent_hmm_list.write("absent_hmms.ids")
    print("%i hmms were not downloaded" % len(absent_hmm_list))

