#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

# TODO: WORKS only for PE DATA and for illumina data
import os
import argparse
from RouToolPa.Tools.Filter import CookiecutterOld
from RouToolPa.Routines.File import check_path, safe_mkdir



parser = argparse.ArgumentParser()

parser.add_argument("-d", "--sample_directory", action="store", dest="samples_dir", required=True,
                    type=lambda s: check_path(os.path.abspath(s)),
                    help="Directory with samples")
parser.add_argument("-s", "--samples", action="store", dest="samples",
                    help="Comma-separated list of subdirectories(one per sample) to handle. "
                         "If not set all subdirectories will be considered as containing samples."
                         "In sample directory should one(in case SE reads) or two(in case PE reads) files."
                         "Filenames should should contain '_1.fq' or '_1.fastq' for forward(left) reads, "
                         " '_2.fq' or '_2.fastq' for reverse(right) reads and '.fq' or '.fastq' for SE reads")
parser.add_argument("-o", "--output_dir", action="store", dest="output_dir",
                    type=lambda s: check_path(os.path.abspath(s)),
                    default="./", help="Directory to write output. Default: current directory")
parser.add_argument("-f", "--fragments_file", action="store", dest="fragments_file", required=True,
                    help="File with k-mers to remove")
parser.add_argument("-c", "--path_to_cookiecutter_dir", action="store", dest="path_to_cookiecutter_dir",
                    default="",
                    help="Path to Cookiecutter directory")

args = parser.parse_args()

samples = args.samples.split(",") if args.samples else sorted(os.listdir(args.samples_dir))
safe_mkdir(args.output_dir)

CookiecutterOld.path = args.path_to_cookiecutter_dir

for sample in samples:
    print("Handling %s" % sample)

    sample_dir = "%s%s/" % (args.samples_dir, sample)

    sample_out_dir = "%s%s/" % (args.output_dir, sample)
    safe_mkdir(sample_out_dir)
    files_from_sample_dir = sorted(os.listdir(sample_dir))

    filtered_files_from_sample_dir = []
    prefix_list = []
    for filename in files_from_sample_dir:
        if ".fq" == filename[-3:]:
            filtered_files_from_sample_dir.append("%s%s" % (sample_dir, filename))
            prefix_list.append("%s%s" % (sample_out_dir, filename[:-3]))
        elif ".fastq" == filename[-6:]:
            filtered_files_from_sample_dir.append("%s%s" % (sample_dir, filename))
            prefix_list.append("%s%s" % (sample_out_dir, filename[:-6]))
        elif ".fq.gz" == filename[-6:]:
            filtered_files_from_sample_dir.append("%s%s" % (sample_dir, filename))
            prefix_list.append("%s%s" % (sample_out_dir, filename[:-6]))
        elif ".fastq.gz" == filename[-9:]:
            filtered_files_from_sample_dir.append("%s%s" % (sample_dir, filename))
            prefix_list.append("%s%s" % (sample_out_dir, filename[:-9]))

    if len(filtered_files_from_sample_dir) % 2 != 0:
        print("Not all read files are paired for sample %s. Skipping..." % sample)
        continue

    number_of_lanes = len(filtered_files_from_sample_dir) / 2

    for lane_number in range(0, number_of_lanes):
        stat_file = "%s.stats" % prefix_list[lane_number*2]
        #output_prefix = "%s%s.TMF" % (sample_out_dir, sample)
        left_reads_file = filtered_files_from_sample_dir[lane_number*2]
        right_reads_file = filtered_files_from_sample_dir[lane_number*2 + 1]

        CookiecutterOld.rm_reads(args.fragments_file, left_reads_file, stat_file, right_reads=right_reads_file,
                                 out_dir=sample_out_dir)