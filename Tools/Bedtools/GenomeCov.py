#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os
from collections import OrderedDict
import numpy as np

from Tools.Abstract import Tool
from Routines import MathRoutines


class GenomeCov(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "bedtools genomecov", path=path, max_threads=max_threads)

    @staticmethod
    def parse_options(input_bam=None, report_zero_coverage=None,
                      one_based_coordinates=None, scale=None, bedgraph_output=False,
                      genome_bed=None):

        options = " -ibam %s" % input_bam
        if bedgraph_output and report_zero_coverage:
            options += " -bga"
        elif bedgraph_output:
            options += " -bg"
        else:
            options += " -d" if one_based_coordinates else " -dz"
        options += " -scale %f" % scale if scale else ""
        options += " -g %s" % genome_bed if genome_bed else ""

        return options

    def get_coverage(self, output, input_bam=None, report_zero_coverage=None,
                     one_based_coordinates=None, scale=None, bedgraph_output=False,
                     genome_bed=None):

        options = self.parse_options(input_bam=input_bam, report_zero_coverage=report_zero_coverage,
                                     one_based_coordinates=one_based_coordinates, scale=scale,
                                     bedgraph_output=bedgraph_output, genome_bed=genome_bed)
        options += " > %s" % output

        self.execute(options)

    def get_coverage_for_gff(self, input_file, genome_bed, output=None):

        options = " -i %s" % input_file
        options += " -g %s" % genome_bed
        options += " > %s" % output if output else ""

        self.execute(options=options)

    def get_bam_coverage_stats(self, input_bam, output_prefix, genome_bed=None, max_coverage=None, min_coverage=None,
                               verbose=True):
        #options_list = []

        each_position_coverage_file = "%s.tab" % output_prefix
        coverage_stat_file = "%s.stat" % output_prefix

        each_position_options = self.parse_options(input_bam=input_bam, report_zero_coverage=True,
                                                   one_based_coordinates=True, genome_bed=genome_bed)
        each_position_options += " > %s" % each_position_coverage_file

        region_options = self.parse_options(input_bam=input_bam, report_zero_coverage=True,
                                            bedgraph_output=True, genome_bed=genome_bed)

        region_options += " > %s.bedgraph" % output_prefix

        options_list = [each_position_options, region_options]

        self.parallel_execute(options_list, threads=2)

        #~/Soft/MAVR/scripts/math/get_stats_from_numer_data.py -i A_ventralis_pe.nodup.q20.tab -o A_ventralis_pe.nodup.q20.stat -l 2

        MathRoutines.get_stats_from_file(each_position_coverage_file, minimum=min_coverage, maximum=max_coverage,
                                         dtype=int, comments="#", delimiter="\t", converters=None, skiprows=0,
                                         usecols=2, unpack=False, ndmin=0, output_file=coverage_stat_file,
                                         verbose=verbose)

    def collapse_coverage_file(self, coverage_file, output_file):

        awk_string = "awk -F'\\t' 'BEGIN {SCAF=\"\"; LEN=\"\"; COV=\"\"} {if (($1 != SCAF)) {if (NR > 1) {printf \"%%s\\t%%s\\t%%s\\n\",SCAF,LEN, COV}; SCAF=$1; LEN=$2; COV=$3} else {LEN=$2; COV=COV\",\"$3}} ; END {printf \"%%s\t%%s\t%%s\n\",SCAF,LEN, COV} %s > %s" % (coverage_file, output_file)

        self.execute(options="", cmd=awk_string)

    @staticmethod
    def extract_data_for_cds_from_collapsed_coverage_file(collapsed_coverage_file, cds_bed_file, output_file,
                                                          skip_transcript_with_no_cds=False):
        cds_dict = OrderedDict()

        with open(cds_bed_file, "r") as cds_fd:
            for line in cds_fd:
                line_list = line.strip().split("\t")
                # convert coordinates to python notation
                cds_dict[line_list[0]] = (int(line_list[1]) - 1, int(line_list[2]))
        with open(collapsed_coverage_file, "r") as col_cov_fd:
            with open(output_file, "w") as out_fd:
                for line in col_cov_fd:
                    transcript_id, length, coverage_array = line.strip().split("\t")

                    if transcript_id not in cds_dict:
                        print("No CDS for transcript %s. %s" % (transcript_id, "Skipping..." if skip_transcript_with_no_cds else ""))
                        if skip_transcript_with_no_cds:
                            continue
                        out_fd.write(line)
                        continue
                    coverage_array = coverage_array.split(",")
                    cds_len = cds_dict[transcript_id][1] - cds_dict[transcript_id][0]
                    out_fd.write("%s\t%s\t%s\n" % (transcript_id, str(cds_len),
                                                   ",".join(coverage_array[cds_dict[transcript_id][0]:cds_dict[transcript_id][1]])))

    @staticmethod
    def analyze_collapsed_coverage_file(collapsed_file, output_file):
        line_number = 0

        def int_through_float(string):
            return int(float(string))

        with open(collapsed_file, "r") as in_fd:
            with open(output_file, "w") as out_fd:

                output_header = "#Scaffold\tLength\tMean_coverage\tMedian_coverage\tMin_coverage\tMax_coverage\tCoveraged_positions\tZero_coverage_position_number\tZero_coveraged_region_number\tLeading_zero_covarage_len\tTrailing_zero_covarage_len\tZero_coverage_region_coordinates\tMean_coverage(without_zerocoveraged_ends)\tMedian_coverage(without_zerocoveraged_ends)\tLength_without_zero_coverage_ends\n"

                out_fd.write(output_header)
                for line in in_fd:
                    line_number += 1
                    #print (line_number)
                    #print [line]
                    if line == "\n" or line == "":    # skip blank lines
                        continue
                    tmp = line.strip().split("\t")
                    #print tmp
                    record_id = tmp[0]
                    record_len = int(tmp[1])

                    try:
                        coverage_array = map(int, tmp[2].split(","))
                    except ValueError:
                        coverage_array = map(int_through_float, tmp[2].split(","))

                    if record_len != len(coverage_array):
                        raise ValueError("Malformed line %i" % line_number)

                    if sum(coverage_array) == 0:
                        continue

                    mean_coverage = float(np.mean(coverage_array))
                    median_coverage = float(np.median(coverage_array))
                    min_coverage = float(np.min(coverage_array))
                    max_coverage = float(np.max(coverage_array))
                    coveraged_position_number = np.count_nonzero(coverage_array)

                    zero_coverage_position_number, zero_coverage_regions_list = MathRoutines.find_flat_regions_in_array(coverage_array, value=0)

                    if zero_coverage_position_number > 0:

                        if zero_coverage_regions_list[0][0] == 0:
                            leading_zero_coverage_len = zero_coverage_regions_list[0][1]
                            start_coverage_coordinate = zero_coverage_regions_list[0][1]
                        else:
                            leading_zero_coverage_len = 0
                            start_coverage_coordinate = 0

                        if zero_coverage_regions_list[-1][0] + zero_coverage_regions_list[-1][1] == record_len:
                            trailing_zero_coverage_len = zero_coverage_regions_list[-1][1]
                            end_coverage_coordinate = zero_coverage_regions_list[-1][0]
                        else:
                            trailing_zero_coverage_len = 0
                            end_coverage_coordinate = record_len

                        zero_coveraged_region_number = len(zero_coverage_regions_list)
                        zero_coverage_coordinates_list = []
                        for (start, length) in zero_coverage_regions_list:
                            zero_coverage_coordinates_list.append("%i-%i" % (start + 1, start + length))

                    else:
                        leading_zero_coverage_len = 0
                        trailing_zero_coverage_len = 0
                        zero_coveraged_region_number = 0
                        start_coverage_coordinate = 0
                        end_coverage_coordinate = record_len
                        zero_coverage_coordinates_list = ["."]

                    length_without_zerocoveraged_ends = end_coverage_coordinate - start_coverage_coordinate
                    mean_coverage_without_zero_coverage_ends = float(np.mean(coverage_array[start_coverage_coordinate:end_coverage_coordinate]))
                    median_coverage_without_zero_coverage_ends = float(np.median(coverage_array[start_coverage_coordinate:end_coverage_coordinate]))

                    out_fd.write("%s\t%i\t%.2f\t%.2f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%s\t%.2f\t%.2f\t%i\n" % (record_id,
                                                                                                           record_len,
                                                                                                           mean_coverage,
                                                                                                           median_coverage,
                                                                                                           min_coverage,
                                                                                                           max_coverage,
                                                                                                           coveraged_position_number,
                                                                                                           zero_coverage_position_number,
                                                                                                           zero_coveraged_region_number,
                                                                                                           leading_zero_coverage_len,
                                                                                                           trailing_zero_coverage_len,
                                                                                                           ",".join(zero_coverage_coordinates_list),
                                                                                                           mean_coverage_without_zero_coverage_ends,
                                                                                                           median_coverage_without_zero_coverage_ends,
                                                                                                           length_without_zerocoveraged_ends
                                                                                                           ))


