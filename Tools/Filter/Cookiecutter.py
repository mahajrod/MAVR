#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os

from Tools.Abstract import Tool


class Cookiecutter(Tool):
    def __init__(self, path="", max_threads=1):
        Tool.__init__(self, "rm_reads", path=path, max_threads=max_threads)

    def filter(self, adapter_file, left_reads, right_reads=None, out_dir="./"):

        options = " --adapters %s" % adapter_file
        options += " -o %s" % out_dir
        options += " -i %s" % left_reads if right_reads is None else ""     # single-end data
        options += " -1 %s -2 %s" % (left_reads, right_reads) if right_reads else ""    # paired-end data

        self.execute(options)

    def extract(self, fragments_file, out_dir, left_reads, right_reads=None):

        options = " --fragments %s" % fragments_file
        options += " -o %s" % out_dir
        options += " -i %s" % left_reads if right_reads is None else ""     # single-end data
        options += " -1 %s -2 %s" % (left_reads, right_reads) if right_reads else ""    # paired-end data

        self.execute(options, cmd="extractor")


if __name__ == "__main__":
    pass