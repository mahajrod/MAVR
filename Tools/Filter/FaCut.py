#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os

from Tools.Abstract import Tool


class FaCut(Tool):
    def __init__(self, path="", max_threads=1):
        Tool.__init__(self, "facut", path=path, max_threads=max_threads)

    def filter_by_mean_quality(self, quality_threshold, forward_reads, reverse_reads,
                               output_prefix, quality_type="phred64", stat_file=None):

        options = " -q %s" % quality_type
        options += " -f %s" % forward_reads
        options += " -r %s" % reverse_reads
        options += " -t %i" % quality_threshold
        options += " -p %s" % output_prefix
        options += " > %s" % stat_file if stat_file else ""

        self.execute(options, cmd="filter_by_mean_quality")


if __name__ == "__main__":
    pass


