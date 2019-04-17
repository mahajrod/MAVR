#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

# TODO: WORKS only for PE DATA and for illumina data
import os
import argparse
from RouToolPa.Tools.Filter import FaCut


#from RouToolPa.Tools.Filter import FastQC

from RouToolPa.Routines import FileRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-d", "--sample_directory", action="store", dest="samples_dir", required=True,
                    type=lambda s: FileRoutines.check_path(os.path.abspath(s)),
                    help="Directory with samples")
parser.add_argument("-s", "--samples", action="store", dest="samples",
                    help="Comma-separated list of subdirectories(one per sample) to handle. "
                         "If not set all subdirectories will be considered as containing samples."
                         "In sample directory should one(in case SE reads) or two(in case PE reads) files."
                         "Filenames should should contain '_1.fq' or '_1.fastq' for forward(left) reads, "
                         " '_2.fq' or '_2.fastq' for reverse(right) reads and '.fq' or '.fastq' for SE reads")
parser.add_argument("-o", "--output_dir", action="store", dest="output_dir",
                    type=lambda s: FileRoutines.check_path(os.path.abspath(s)),
                    default="./", help="Directory to write output. Default: current directory")
#parser.add_argument("-t", "--threads", action="store", dest="threads", default=1, type=int,
#                    help="Number of threads to use in Trimmomatic. Default - 1.")
parser.add_argument("-q", "--average_quality_threshold", action="store", dest="average_quality_threshold", default=15,
                    type=int,
                    help="Mean quality threshold for read. Works only if -q/--average_quality_threshold is set"
                         "Default - 15.")
parser.add_argument("-u", "--score_type", action="store", dest="score_type", default="phred64",
                    help="Phred quality score type. Allowed: phred33, phred64. Default: phred64")
parser.add_argument("-n", "--name_type", action="store", dest="name_type", default="short",
                    help="Type of read name. Required to gather per tile filtering statistics. Default: short")
args = parser.parse_args()

samples = args.samples.split(",") if args.samples else sorted(os.listdir(args.samples_dir))
FileRoutines.safe_mkdir(args.output_dir)
for sample in samples:
    print("Handling %s" % sample)

    sample_dir = "%s%s/" % (args.samples_dir, sample)

    sample_out_dir = "%s%s/" % (args.output_dir, sample)
    FileRoutines.safe_mkdir(sample_out_dir)
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
        #output_prefix = "%s%s.TMF" % (sample_out_dir, sample)
        left_reads_file = filtered_files_from_sample_dir[lane_number*2]
        right_reads_file = filtered_files_from_sample_dir[lane_number*2 + 1]

        """
        filter_string = "time filter_by_mean_quality -t %i -f %s -r %s -q %s -p %s" % (args.average_quality_threshold,
                                                                                       left_reads_file,
                                                                                       right_reads_file,
                                                                                       args.score_type,
                                                                                       prefix_list[lane_number*2])
        """
        #print(left_reads_file)
        #print(right_reads_file)
        #print(args.score_type)
        #print(args.average_quality_threshold)
        #print(prefix_list[lane_number*2])

        #print filter_string
        #os.system(filter_string)
        FaCut.timelog = "%s.timelog" % prefix_list[lane_number*2]
        FaCut.filter_by_mean_quality(args.average_quality_threshold, left_reads_file, right_reads_file,
                                     prefix_list[lane_number*2], quality_type=args.score_type,
                                     stat_file="%s.stat" % prefix_list[lane_number*2],
                                     name_type=args.name_type)

