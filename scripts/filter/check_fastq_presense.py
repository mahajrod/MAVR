#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

# TODO: WORKS only for PE DATA and for illumina data
import os
import sys
import argparse

from Tools.Filter import FaCut
#from Tools.Filter import FastQC

from Routines.File import check_path, save_mkdir

def check_if_files_form_pair(first_filename, second_filename, name_type):
    first_list = first_filename.split("_")
    second_list = second_filename("_")
    if len(first_list) != len(second_list):
        return False

    if name_type == "botswana":
        for j in range(0, len(first_list) - 1):
            if first_list[j] != second_list[j]:
                return False
        if first_list[-1][1:] != second_list[-1][1:]:
            return False
        if first_list[-1][1] != 1 or second_list[-1][1] != 2:
            return False
        return True


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
parser.add_argument("-o", "--output", action="store", dest="output",
                    type=lambda s: check_path(os.path.abspath(s)), requireed=True,
                    help="File to write output.")
#parser.add_argument("-t", "--threads", action="store", dest="threads", default=1, type=int,
#                    help="Number of threads to use in Trimmomatic. Default - 1.")

parser.add_argument("-n", "--name_type", action="store", dest="name_type", default="botswana",
                    help="Type of file name. Allowed: botswana"
                         "Botswana 121212_I595_FCC1H08ACXX_L8_SZAIPI018897-93_1.fq. Default: botswana")
args = parser.parse_args()

samples = args.samples.split(",") if args.samples else sorted(os.listdir(args.samples_dir))
save_mkdir(args.output_dir)

out_fd = open(arg)
for sample in samples:
    print("Handling %s" % sample)

    sample_dir = "%s%s/" % (args.samples_dir, sample)

    sample_out_dir = "%s%s/" % (args.output_dir, sample)
    save_mkdir(sample_out_dir)
    files_from_sample_dir = sorted(os.listdir(sample_dir))

    read_files_list = []
    splited_read_files = []

    for filename in files_from_sample_dir:
        if ".fq" == filename[-3:] or ".fastq" == filename[-6:]  \
                or ".fq.gz" == filename[-6:] or ".fastq.gz" == filename[-9:]:
            read_files_list.append(filename)

    number_of_files = len(read_files_list)
    i = 0
    number_of_paires_of_files = 0
    number_of_unpaired_files = 0
    while i < number_of_files:
        if i+1 != number_of_files:
            if check_if_files_form_pair(read_files_list[i], read_files_list[i+1]):
                number_of_paires_of_files += 1
                i += 2
            else:
                number_of_unpaired_files += 1
                i += 1
        else:
            number_of_unpaired_files += 1
            i += 1

    if len(read_files_list) % 2 != 0:
        print("Not all read files are paired for sample %s. Skipping..." % sample)
        continue



    for lane_number in range(0, number_of_lanes):
        #output_prefix = "%s%s.TMF" % (sample_out_dir, sample)
        left_reads_file = filtered_files_from_sample_dir[lane_number*2]
        right_reads_file = filtered_files_from_sample_dir[lane_number*2 + 1]

        """
        filter_string = "time filter_by_mean_quality -t %i -f %s -r %s -q %s -p %s" % (args.average_quality_threshold,
                                                                                       left_reads_file,
                                                                                       right_reads_file,
                                                                                       args.score_type,
                                                                                       read_files_list[lane_number*2])
        """
        #print(left_reads_file)
        #print(right_reads_file)
        #print(args.score_type)
        #print(args.average_quality_threshold)
        #print(read_files_list[lane_number*2])

        #print filter_string
        #os.system(filter_string)
        FaCut.timelog = "%s.timelog" % read_files_list[lane_number*2]
        FaCut.filter_by_mean_quality(args.average_quality_threshold, left_reads_file, right_reads_file,
                                     read_files_list[lane_number*2], quality_type=args.score_type,
                                     stat_file="%s.stat" % read_files_list[lane_number*2],
                                     name_type=args.name_type)

