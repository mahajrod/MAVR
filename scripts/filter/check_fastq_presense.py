#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os
import argparse
from RouToolPa.Routines.File import check_path


def check_if_files_form_pair(first_filename, second_filename, name_type):
    first_list = first_filename.split("_")
    second_list = second_filename.split("_")
    if len(first_list) != len(second_list):
        return False

    if name_type == "botswana":
        for j in range(0, len(first_list) - 1):
            if first_list[j] != second_list[j]:
                return False
        if first_list[-1][1:] != second_list[-1][1:]:
            return False
        if first_list[-1][0] != "1" or second_list[-1][0] != "2":
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
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="File to write output.")
parser.add_argument("-n", "--name_type", action="store", dest="name_type", default="botswana",
                    help="Type of file name. Default: botswana")

args = parser.parse_args()

samples = args.samples.split(",") if args.samples else sorted(os.listdir(args.samples_dir))


out_fd = open(args.output, "w")
out_fd.write("#Sample\tNumber_of_files\tNumber_of_paires\tNumber_of_unpaired_files\n")

for sample in samples:
    print("Handling %s" % sample)

    sample_dir = "%s%s/" % (args.samples_dir, sample)

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
            if check_if_files_form_pair(read_files_list[i], read_files_list[i+1], args.name_type):
                number_of_paires_of_files += 1
                i += 2
            else:
                number_of_unpaired_files += 1
                i += 1
        else:
            number_of_unpaired_files += 1
            i += 1

    number_of_files = len(read_files_list)
    if number_of_files % 2 != 0:
        print("Not all read files are paired for sample %s. Skipping..." % sample)

    out_fd.write("%s\t%i\t%i\t%i\n" % (sample, number_of_files,
                                       number_of_paires_of_files, number_of_unpaired_files))

