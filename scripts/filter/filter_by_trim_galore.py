#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse
from RouToolPa.Tools.Filter import TrimGalore
from RouToolPa.Tools.Filter import FastQC

parser = argparse.ArgumentParser()

parser.add_argument("-s", "--sample_directory", action="store", dest="sample_dir",
                    help="Directory with samples")
parser.add_argument("-d", "--description_file", action="store", dest="description",
                    help="File with description of samples")
parser.add_argument("-t", "--trim_subdir", action="store", dest="trim_subdir", default="trimmed",
                    help="Subdirectory to write filtered data")
parser.add_argument("-l", "--min_length", action="store", dest="min_len", default=30, type=int,
                    help="Minimum length of read to retain")
parser.add_argument("-q", "--quality_threshold", action="store", dest="quality_threshold", default=20, type=int,
                    help="Quality threshold of reads")
args = parser.parse_args()

with open(args.description, "r") as desr_fd:
    header = desr_fd.readline().strip().split("\t")
    for line in desr_fd:
        sample_name, platform, library, library_type, quality_type, read_length, number_of_pairs, left_reads_trim, right_reads_trim, \
        seq_to_trim, total_bases, coverage, note = line.strip().split("\t")
        os.chdir(args.sample_dir)
        os.chdir(sample_name)
        if library_type == "PE":
            left_reads = "%s_1.fastq" % sample_name
            right_reads = "%s_2.fastq" % sample_name
        else:
            left_reads = "%s.fastq" % sample_name
            right_reads = None
        left_reads_trim = None if left_reads_trim == "." else int(left_reads_trim)
        right_reads_trim = None if right_reads_trim == "." else int(right_reads_trim)
        TrimGalore.filter(args.min_len, left_reads, left_reads_trim,
                          reverse_reads=right_reads, reverse_trim=right_reads_trim,
                          quality_score=quality_type, adapter=seq_to_trim,
                          quality_treshold=args.quality_threshold, output_folder=args.trim_subdir)
        os.chdir(args.trim_subdir)
        FastQC.threads = 2
        FastQC.check("*.f*q")

