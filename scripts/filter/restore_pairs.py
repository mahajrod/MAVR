#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
"""
Script is very slow because it uses Biopython SeqIO to parse fastq.
Use it only for small fastqs
"""
import os
import re
import sys
import argparse

from Bio import SeqIO

from RouToolPa.Routines import SequenceRoutines


def get_list_from_string(s):
    return s.split(",")

parser = argparse.ArgumentParser()

parser.add_argument("-l", "--input_left", action="store", dest="input_left", type=get_list_from_string, required=True,
                    help="Comma-separated list of files with left reads")
parser.add_argument("-r", "--input_right", action="store", dest="input_right", type=get_list_from_string, required=True,
                    help="Comma-separated list of files with left reads")
parser.add_argument("-o", "--out_prefix", action="store", dest="out_prefix", required=True,
                    help="Prefix of output files")

args = parser.parse_args()

out_left = "%s_1.fastq" % args.out_prefix
out_right = "%s_2.fastq" % args.out_prefix
out_left_se = "%s_1.se.fastq" % args.out_prefix
out_right_se = "%s_2.se.fastq" % args.out_prefix

out_left_fd = open("%s_1.fastq" % args.out_prefix, "w")
out_right_fd = open("%s_2.fastq" % args.out_prefix, "w")
out_left_se_fd = open("%s_1.se.fastq" % args.out_prefix, "w")
out_right_se_fd = open("%s_2.se.fastq" % args.out_prefix, "w")

left_input_reads_dict = SeqIO.index_db("left_in_reads.idx", args.input_left, "fastq")
right_input_reads_dict = SeqIO.index_db("right_in_reads.idx", args.input_right, "fastq")

left_input_set = set(left_input_reads_dict.keys())
right_input_set = set(right_input_reads_dict.keys())


SeqIO.write(SequenceRoutines.record_by_id_generator(left_input_reads_dict,
                                                    sorted(left_input_set & right_input_set),
                                                    verbose=True), out_left, "fastq")
SeqIO.write(SequenceRoutines.record_by_id_generator(right_input_reads_dict,
                                                    sorted(left_input_set & right_input_set),
                                                    verbose=True), out_right, "fastq")
SeqIO.write(SequenceRoutines.record_by_id_generator(left_input_reads_dict,
                                                    left_input_set - right_input_set,
                                                    verbose=True), out_left_se, "fastq")
SeqIO.write(SequenceRoutines.record_by_id_generator(right_input_reads_dict,
                                                    right_input_set - left_input_set,
                                                    verbose=True), out_right_se, "fastq")
out_left_fd.close()
out_right_fd.close()
out_left_se_fd.close()
out_right_se_fd.close()
os.remove("left_in_reads.idx")
os.remove("right_in_reads.idx")
