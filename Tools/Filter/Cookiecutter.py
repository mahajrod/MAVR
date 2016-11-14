#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os

from Tools.Abstract import Tool


class Cookiecutter(Tool):
    def __init__(self, path="", max_threads=1):
        Tool.__init__(self, "cookiecutter", path=path, max_threads=max_threads)

    def parse_common_options(self, adapter_file, left_reads, right_reads=None, out_dir="./"):

        options = " --fragments %s" % adapter_file
        #options += " -t %i" % self.threads
        options += " -o %s" % out_dir
        options += " -i %s" % left_reads if right_reads is None else ""     # single-end data
        options += " -1 %s -2 %s" % (left_reads, right_reads) if right_reads else ""    # paired-end data

        return options

    def remove(self, adapter_file, left_reads, right_reads=None, out_dir="./"):
        options = self.parse_common_options(adapter_file, left_reads, right_reads=right_reads, out_dir=out_dir)
        self.execute(options, cmd="remove")

    def extract(self, adapter_file, left_reads, right_reads=None, out_dir="./"):
        options = self.parse_common_options(adapter_file, left_reads, right_reads=right_reads, out_dir=out_dir)
        self.execute(options, cmd="extract")

    def separate(self, adapter_file, left_reads, right_reads=None, out_dir="./"):
        options = self.parse_common_options(adapter_file, left_reads, right_reads=right_reads, out_dir=out_dir)
        self.execute(options, cmd="separate")

    def rm_reads(self, adapter_file, left_reads, stats_file, right_reads=None, out_dir="./", use_dust_filter=False,
                 dust_cutoff=None, dust_window_size=None, use_N_filter=False,
                 read_length_cutoff=None, polyGC_length_cutoff=None):
        options = self.parse_common_options(adapter_file, left_reads, right_reads=right_reads, out_dir=out_dir)
        options += " -p %i" % polyGC_length_cutoff if polyGC_length_cutoff else ""
        options += " -l %i" % read_length_cutoff if read_length_cutoff else ""
        options += " -d" if use_dust_filter else ""
        options += " -c %i" % dust_cutoff if dust_cutoff else ""
        options += " -k %i" % dust_window_size if dust_window_size else ""
        options += " -N" if use_N_filter else ""
        options += " > %s" % stats_file

        self.execute(options, cmd="rm_reads")


class CookiecutterOld(Tool):
    def __init__(self, path="", max_threads=1):
        Tool.__init__(self, "cookiecutter", path=path, max_threads=max_threads)

    def rm_reads(self, adapter_file, left_reads, stats_file, right_reads=None, out_dir="./"):
        options = " --fragments %s" % adapter_file
        options += " -o %s" % out_dir
        options += " -i %s" % left_reads if right_reads is None else ""     # single-end data
        options += " -1 %s -2 %s" % (left_reads, right_reads) if right_reads else ""    # paired-end data
        options += " > %s" % stats_file

        self.execute(options, cmd="rm_reads")

if __name__ == "__main__":
    pass


