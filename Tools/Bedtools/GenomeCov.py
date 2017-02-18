#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os

from Tools.Abstract import Tool


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

    def convert_coverage_file(self, coverage_file, output_file):
        pass


if __name__ == "__main__":
    pass
