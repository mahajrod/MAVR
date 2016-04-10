#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os

from Tools.Abstract import Tool


class DiscovarDeNovo(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "DiscovarDeNovo", path=path, max_threads=max_threads)

    def get_stat(self, list_of_read_files, output_dir, reference):

        options = " NUM_THREADS=%i" % self.threads
        options += " READS=%s" % ",".join(list_of_read_files)
        options += " OUT_DIR=%s" % output_dir
        options += " REFHEAD=%s" % reference

        self.execute(options)

if __name__ == "__main__":
    pass
