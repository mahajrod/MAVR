#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os
import numpy as np

from Tools.Abstract import Tool
from Routines import MathRoutines

class GenomeCov(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "bedtools genomecov", path=path, max_threads=max_threads)

    def get_coverage(self, input_bam, genome_bed, output, report_zero_coverage=None, one_based_coordinates=True, scale=None):

        options = " -ibam %s" % input_bam
        # options += " -bga" if report_zero_coverage else " -bg"
        # options += " -bga" if report_zero_coverage else " -bg"
        options += " -d" if one_based_coordinates else " -dz"
        options += " -scale %f" % scale if scale else ""
        options += " -g %s" % genome_bed
        options += " > %s" % output

        self.execute(options)

    def collapse_coverage_file(self, coverage_file, output_file):

        awk_string = "awk -F'\\t' 'BEGIN {SCAF=\"\"; LEN=\"\"; COV=\"\"} {if ($1 != SCAF) {printf \"%s\\t%s\\t%s\\n\",SCAF,LEN,COV; SCAF=$1; LEN=$2; COV=$3} else {LEN=$2; COV=COV\",\"$3}} ; END {printf \"%s\\t%s\\t%s\\n\",SCAF,LEN,COV}' %s > %s" % (coverage_file, output_file)

        self.execute(options="", cmd=awk_string)

    @staticmethod
    def analyze_collapsed_coverage_file(collapsed_file, output_file):
        line_number = 0
        with open(collapsed_file, "r") as in_fd:
            with open(output_file, "w") as out_fd:

                output_header = "#Scaffold\tLength\tMean_coverage\tMedian_coverage\tMin_coverage\tMax_coverage\tCoveraged_positions\tZero_coverage_position_number\tZero_coveraged_region_number\tLeading_zero_covarage_len\tTrailing_zero_covarage_len\tZero_coverage_region_coordinates\n"

                out_fd.write(output_header)
                for line in in_fd:
                    line_number += 1
                    print (line_number)
                    print [line]
                    if line == "\n" or line == "":    # skip blank lines
                        continue
                    tmp = line.strip().split("\t")
                    #print tmp
                    record_id = tmp[0]
                    record_len = int(tmp[1])
                    coverage_array = map(int, tmp[2].split(","))

                    if record_len != len(coverage_array):
                        raise ValueError("Malformed line %i" % line_number)

                    mean_covarage = float(np.mean(coverage_array))
                    median_coverage = float(np.median(coverage_array))
                    min_coverage = float(np.min(coverage_array))
                    max_coverage = float(np.max(coverage_array))
                    coveraged_position_number = np.count_nonzero(coverage_array)

                    zero_coverage_position_number, zero_coverage_regions_list = MathRoutines.find_flat_regions_in_array(coverage_array, value=0)

                    if zero_coverage_position_number > 0:
                        leading_zero_covarage_len = zero_coverage_regions_list[0][1]
                        trailing_zero_covarage_len = zero_coverage_regions_list[-1][1]
                        zero_coveraged_region_number = len(zero_coverage_regions_list)
                        zero_coverage_coordinates_list = []
                        for (start, length) in zero_coverage_regions_list:
                            zero_coverage_coordinates_list.append("%i-%i" % (start + 1, start + length))

                    else:
                        leading_zero_covarage_len = 0
                        trailing_zero_covarage_len = 0
                        zero_coveraged_region_number = 0
                        zero_coverage_coordinates_list = ["."]

                    out_fd.write("%s\t%i\t%.2f\t%.2f\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%s\n" % (record_id,
                                                                                           record_len,
                                                                                           mean_covarage,
                                                                                           median_coverage,
                                                                                           min_coverage,
                                                                                           max_coverage,
                                                                                           coveraged_position_number,
                                                                                           zero_coverage_position_number,
                                                                                           zero_coveraged_region_number,
                                                                                           leading_zero_covarage_len,
                                                                                           trailing_zero_covarage_len,
                                                                                           ",".join(zero_coverage_coordinates_list)))


