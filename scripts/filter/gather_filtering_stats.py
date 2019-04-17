#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os
import argparse
import numpy as np
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
"""
#parser.add_argument("-t", "--threads", action="store", dest="threads", default=1, type=int,
#                    help="Number of threads to use in Trimmomatic. Default - 1.")
parser.add_argument("-q", "--average_quality_threshold", action="store", dest="average_quality_threshold", default=15,
                    type=int,
                    help="Quality threshold for sliding window. Works only if -q/--average_quality_threshold is set"
                         "Default - 15.")
parser.add_argument("-u", "--score_type", action="store", dest="score_type", default="phred64",
                    help="Phred quality score type. Allowed: phred33, phred64. Default: phred64")
parser.add_argument("-n", "--name_type", action="store", dest="name_type", default="short",
                    help="Type of read name. Required to gather per tile filtering statistics. Default: short")
"""
args = parser.parse_args()

samples = args.samples.split(",") if args.samples else sorted(os.listdir(args.samples_dir))
FileRoutines.safe_mkdir(args.output_dir)

overall_stat_file = "%s/overall_samples.stat" % args.output_dir
overall_stat_fd = open(overall_stat_file, "w")
overall_stat_fd.write("#Sample_id\tTotal_pairs\tRetained_pairs\tRetained_pairs_percent\tMin_pairs_retained_in_tiles\n")

for sample in samples:
    print("Handling %s" % sample)

    sample_dir = "%s%s/" % (args.samples_dir, sample)

    sample_out_dir = "%s%s/" % (args.output_dir, sample)
    FileRoutines.safe_mkdir(sample_out_dir)
    files_from_sample_dir = sorted(os.listdir(sample_dir))

    stat_files_from_sample_dir = []
    prefix_list = []
    for filename in files_from_sample_dir:
        if ".stat" == filename[-5:]:
            stat_files_from_sample_dir.append("%s%s" % (sample_dir, filename))
            prefix_list.append("%s%s" % (sample_out_dir, filename[:-5]))

    number_of_stat_files = len(stat_files_from_sample_dir)

    percent_total_sample_stat_file = "%s/%s.sample.percent.stats" % (sample_out_dir, sample)

    sample_stat_fd = open(percent_total_sample_stat_file, "w")
    sample_stat_fd.write("Read_group\t"
                         "Paires_retained\tForward_only_retained\tReverse_only_retained\tPairs_discarded\t"
                         "Tile_min_pairs_retained\tTile_min_forward_only_retained\t"
                         "Tile_min_reverse_only_retained\tTile_min_pairs_discarded\t"
                         "Tile_max_pairs_retained\tTile_max_forward_only_retained\t"
                         "Tile_max_reverse_only_retained\tTile_max_pairs_discarded\t"
                         "Tile_mean_pairs_retained\tTile_mean_forward_only_retained\t"
                         "Tile_mean_reverse_only_retained\tTile_mean_pairs_discarded\t"
                         "Tile_median_pairs_retained\tTile_median_forward_only_retained\t"
                         "Tile_median_reverse_only_retained\tTile_median_pairs_discarded\n")

    total_reads_sample_stats = []
    min_percent_retained_pairs_in_tile_list = []
    for stat_file_index in range(0, number_of_stat_files):

        percent_total_stat_file = "%s.total.percent.stats" % prefix_list[stat_file_index]
        percent_tile_stat_file = "%s.tile.percent.stats" % prefix_list[stat_file_index]

        tile_description_list = []
        total_reads_list = []
        tile_stats_list = []

        with open(stat_files_from_sample_dir[stat_file_index], "r") as stat_fd:
            try:
                line = stat_fd.next()

                while line[:13] != "instrument_id":
                    if line[:15] == "Paires retained" or line[:14] == "Pairs retained":
                        pairs_retained = float(line.strip().split("\t")[-1])
                    elif line[:21] == "Forward only retained":
                        forward_only_retained = float(line.strip().split("\t")[-1])
                    elif line[:21] == "Reverse only retained":
                        reverse_only_retained = float(line.strip().split("\t")[-1])
                    elif line[:15] == "Pairs discarded":
                        pairs_discarded = float(line.strip().split("\t")[-1])
                    line = stat_fd.next()
                line = stat_fd.next()
                total_stats_list = np.array([pairs_retained, forward_only_retained, reverse_only_retained, pairs_discarded])
                for line in stat_fd:
                    line_list = line.strip().split("\t")
                    tile_stats = map(float, line_list[-4:])
                    # skip absent tile
                    if sum(tile_stats) == 0:
                        print("\tTile %s is absent in input data for %s" % (line_list[4], prefix_list[stat_file_index]))
                        continue
                    tile_description_list.append(line_list[:-4])
                    tile_stats_list.append(tile_stats)
                    total_reads_list.append(sum(tile_stats))

            except StopIteration:
                print("\tEmpty .stat file for %s" % prefix_list[stat_file_index])
                continue

        total_reads_sample_stats.append(total_stats_list)
        total_reads_list = np.array(total_reads_list)
        # tile_stats
        tile_stats_list = np.array(tile_stats_list)
        percent_stats_list = tile_stats_list / total_reads_list[:, None]

        #print(percent_stats_list)

        # total_stats

        total_percent_stats = total_stats_list / sum(total_stats_list)
        samples_mean_percent_stats = np.mean(percent_stats_list, axis=0)
        samples_median_percent_stats = np.median(percent_stats_list, axis=0)
        samples_max_percent_stats = np.max(percent_stats_list, axis=0)
        samples_min_percent_stats = np.min(percent_stats_list, axis=0)

        min_percent_retained_pairs_in_tile_list.append(samples_min_percent_stats[0])

        with open(percent_total_stat_file, "w") as percent_stats_fd:
            percent_stats_fd.write("Paires retained\tForward only retained\tReverse only retained\tPairs discarded\n")
            percent_stats_fd.write("%s\n" % "\t".join(map(lambda f: "%.3f" % f, total_percent_stats)))

        sample_stat_fd.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (prefix_list[stat_file_index].split("/")[-1],
                                                           "\t".join(map(lambda f: "%.3f" % f, total_percent_stats)),
                                                           "\t".join(map(lambda f: "%.3f" % f, samples_min_percent_stats)),
                                                           "\t".join(map(lambda f: "%.3f" % f, samples_max_percent_stats)),
                                                           "\t".join(map(lambda f: "%.3f" % f, samples_mean_percent_stats)),
                                                           "\t".join(map(lambda f: "%.3f" % f, samples_median_percent_stats))))

    total_reads_sample_stats = np.array(total_reads_sample_stats)
    # print(total_reads_sample_stats)

    summed_reads_sample_stats = np.sum(total_reads_sample_stats, axis=0)
    total_number_of_pairs = np.sum(summed_reads_sample_stats)
    percent_summed_reads_sample_stats = summed_reads_sample_stats/total_number_of_pairs
    sample_stat_fd.write("Total\t%s\n" % ("\t".join(map(lambda f: "%.3f" % f, percent_summed_reads_sample_stats))))
    sample_stat_fd.close()

    overall_stat_fd.write("%s\t%.0f\t%.0f\t%.3f\t%.3f\n" % (sample, total_number_of_pairs,
                                                            summed_reads_sample_stats[0],
                                                            percent_summed_reads_sample_stats[0],
                                                            min(min_percent_retained_pairs_in_tile_list)))


overall_stat_fd.close()








